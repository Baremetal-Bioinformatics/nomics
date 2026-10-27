use super::cigar::dangling_cigar;
use crate::jmath;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet, VecDeque};
use std::hash::{BuildHasherDefault, Hasher};

#[derive(Default)]
struct FxHasher(u64);

impl FxHasher {
    fn add(&mut self, value: u64) {
        self.0 = (self.0.rotate_left(5) ^ value).wrapping_mul(0x517c_c1b7_2722_0a95);
    }
}

impl Hasher for FxHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, mut bytes: &[u8]) {
        while let Some((chunk, rest)) = bytes.split_first_chunk::<8>() {
            self.add(u64::from_ne_bytes(*chunk));
            bytes = rest;
        }
        let mut tail = [0_u8; 8];
        tail[..bytes.len()].copy_from_slice(bytes);
        self.add(u64::from_ne_bytes(tail));
    }

    fn write_u64(&mut self, value: u64) {
        self.add(value);
    }

    fn write_usize(&mut self, value: usize) {
        self.add(value as u64);
    }
}

// These tables are membership-only; vertex and edge vectors remain the sole
// source of traversal order and first-edge tie-breaking.
type FastMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;
type FastSet<K> = HashSet<K, BuildHasherDefault<FxHasher>>;

fn base_bits(base: u8) -> Option<u64> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn packed_kmer(sequence: &[u8]) -> Option<u64> {
    if sequence.len() > 32 {
        return None;
    }
    sequence.iter().try_fold(0_u64, |packed, base| {
        let bits = base_bits(*base)?;
        Some((packed << 2) | bits)
    })
}

struct KmerSet {
    kmer: usize,
    packed: FastSet<u64>,
    other: FastSet<Vec<u8>>,
}

impl KmerSet {
    fn new(kmer: usize) -> Self {
        Self {
            kmer,
            packed: FastSet::default(),
            other: FastSet::default(),
        }
    }

    fn contains(&self, sequence: &[u8]) -> bool {
        if sequence.len() == self.kmer
            && let Some(packed) = packed_kmer(sequence)
        {
            self.packed.contains(&packed)
        } else {
            self.other.contains(sequence)
        }
    }

    fn insert(&mut self, sequence: &[u8]) -> bool {
        if sequence.len() == self.kmer
            && let Some(packed) = packed_kmer(sequence)
        {
            self.packed.insert(packed)
        } else {
            self.other.insert(sequence.to_vec())
        }
    }

    fn with_capacity(kmer: usize, capacity: usize) -> Self {
        Self {
            kmer,
            packed: FastSet::with_capacity_and_hasher(capacity, Default::default()),
            other: FastSet::default(),
        }
    }

    fn clear(&mut self) {
        self.packed.clear();
        self.other.clear();
    }

    fn insert_packed(&mut self, packed: u64) -> bool {
        self.packed.insert(packed)
    }

    fn len(&self) -> usize {
        self.packed.len() + self.other.len()
    }
}

// Graph construction only queries fixed-length kmers. Canonical kmers up to
// 32 bases therefore have an exact, allocation-free key; all others retain
// bytewise equality. Neither table is iterated, so hashing cannot affect order.
struct KmerMap {
    kmer: usize,
    packed: FastMap<u64, VertexId>,
    other: FastMap<Vec<u8>, VertexId>,
}

impl KmerMap {
    fn with_capacity(kmer: usize, capacity: usize) -> Self {
        Self {
            kmer,
            packed: FastMap::with_capacity_and_hasher(capacity, Default::default()),
            other: FastMap::default(),
        }
    }

    fn get(&self, sequence: &[u8]) -> Option<&VertexId> {
        if sequence.len() == self.kmer
            && let Some(packed) = packed_kmer(sequence)
        {
            self.packed.get(&packed)
        } else {
            self.other.get(sequence)
        }
    }

    fn remove(&mut self, sequence: &[u8]) {
        if sequence.len() == self.kmer
            && let Some(packed) = packed_kmer(sequence)
        {
            self.packed.remove(&packed);
        } else {
            self.other.remove(sequence);
        }
    }

    fn insert_if_absent(&mut self, sequence: &[u8], vertex: VertexId) {
        if sequence.len() == self.kmer
            && let Some(packed) = packed_kmer(sequence)
        {
            self.packed.entry(packed).or_insert(vertex);
        } else {
            self.other.entry(sequence.to_vec()).or_insert(vertex);
        }
    }

    fn len(&self) -> usize {
        self.packed.len() + self.other.len()
    }
}

type VertexId = usize;
type EdgeId = usize;

#[derive(Clone)]
struct Vertex {
    sequence: Vec<u8>,
    /// Last base of `sequence` at creation; read only while threading,
    /// before any sequence is rewritten.
    last: u8,
    alive: bool,
    incoming: Vec<EdgeId>,
    outgoing: Vec<EdgeId>,
    live_incoming: usize,
    live_outgoing: usize,
}

#[derive(Clone)]
struct Edge {
    source: VertexId,
    target: VertexId,
    multiplicity: i32,
    current_sample: i32,
    pruning_multiplicity: i32,
    is_ref: bool,
    alive: bool,
}

struct Graph {
    kmer: usize,
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
    unique: KmerMap,
    non_unique: KmerSet,
    ref_source_kmer: Vec<u8>,
    ref_source: VertexId,
    ref_sink: VertexId,
    live_vertices: usize,
    /// kmerToVertexMap is only read while threading and by isLowQualityGraph;
    /// after that check its maintenance is unobservable and is skipped.
    track_unique: bool,
}

impl Graph {
    fn new(kmer: usize, non_unique: KmerSet, capacity: usize) -> Self {
        Self {
            kmer,
            vertices: Vec::with_capacity(capacity),
            edges: Vec::with_capacity(capacity),
            unique: KmerMap::with_capacity(kmer, capacity),
            non_unique,
            ref_source_kmer: Vec::new(),
            ref_source: usize::MAX,
            ref_sink: usize::MAX,
            live_vertices: 0,
            track_unique: true,
        }
    }

    fn add_vertex(&mut self, sequence: Vec<u8>) -> VertexId {
        let id = self.vertices.len();
        if self.track_unique && !self.non_unique.contains(&sequence) {
            self.unique.insert_if_absent(&sequence, id);
        }
        let last = sequence.last().copied().unwrap_or(0);
        self.vertices.push(Vertex {
            last,
            sequence,
            alive: true,
            incoming: Vec::new(),
            outgoing: Vec::new(),
            live_incoming: 0,
            live_outgoing: 0,
        });
        self.live_vertices += 1;
        id
    }

    fn vertex_for(&mut self, kmer: &[u8], allow_ref_source: bool) -> VertexId {
        if (allow_ref_source || kmer != self.ref_source_kmer)
            && let Some(id) = self.unique.get(kmer)
        {
            return *id;
        }
        self.add_vertex(kmer.to_vec())
    }

    fn live_outgoing(&self, vertex: VertexId) -> Vec<EdgeId> {
        self.vertices[vertex]
            .outgoing
            .iter()
            .copied()
            .filter(|edge| self.edges[*edge].alive)
            .collect()
    }

    fn live_outgoing_iter(&self, vertex: VertexId) -> impl Iterator<Item = EdgeId> + '_ {
        self.vertices[vertex]
            .outgoing
            .iter()
            .copied()
            .filter(|edge| self.edges[*edge].alive)
    }

    fn live_incoming(&self, vertex: VertexId) -> Vec<EdgeId> {
        self.vertices[vertex]
            .incoming
            .iter()
            .copied()
            .filter(|edge| self.edges[*edge].alive)
            .collect()
    }

    fn live_incoming_iter(&self, vertex: VertexId) -> impl Iterator<Item = EdgeId> + '_ {
        self.vertices[vertex]
            .incoming
            .iter()
            .copied()
            .filter(|edge| self.edges[*edge].alive)
    }

    fn add_edge(&mut self, source: VertexId, target: VertexId, is_ref: bool, count: i32) -> EdgeId {
        if let Some(existing) = self
            .live_outgoing_iter(source)
            .find(|edge| self.edges[*edge].target == target)
        {
            return existing;
        }
        let id = self.edges.len();
        self.edges.push(Edge {
            source,
            target,
            multiplicity: count,
            current_sample: count,
            pruning_multiplicity: count,
            is_ref,
            alive: true,
        });
        self.vertices[source].outgoing.push(id);
        self.vertices[target].incoming.push(id);
        self.vertices[source].live_outgoing += 1;
        self.vertices[target].live_incoming += 1;
        id
    }

    fn add_base_edge(
        &mut self,
        source: VertexId,
        target: VertexId,
        is_ref: bool,
        count: i32,
    ) -> EdgeId {
        // GATK's graphs are jgrapht DefaultDirectedGraphs: addEdge between a
        // pair that is already connected is refused, and the existing edge
        // (with its multiplicity) stands. SharedSequenceMerger hits this
        // when two merged predecessors share a source.
        if let Some(existing) = self
            .live_outgoing_iter(source)
            .find(|edge| self.edges[*edge].target == target)
        {
            return existing;
        }
        let id = self.edges.len();
        self.edges.push(Edge {
            source,
            target,
            multiplicity: count,
            current_sample: 0,
            pruning_multiplicity: count,
            is_ref,
            alive: true,
        });
        self.vertices[source].outgoing.push(id);
        self.vertices[target].incoming.push(id);
        self.vertices[source].live_outgoing += 1;
        self.vertices[target].live_incoming += 1;
        id
    }

    fn flush_sample(&mut self) {
        for edge in &mut self.edges {
            if edge.alive {
                edge.pruning_multiplicity = edge.pruning_multiplicity.max(edge.current_sample);
                edge.current_sample = 0;
            }
        }
    }

    fn increment_backwards(&mut self, sequence: &[u8], vertex: VertexId, offset: isize) {
        if offset < 0 {
            return;
        }
        let mut incoming = self.live_incoming_iter(vertex);
        let Some(edge) = incoming.next() else { return };
        if incoming.next().is_some() {
            return;
        }
        drop(incoming);
        let previous = self.edges[edge].source;
        if self.vertices[previous].last == sequence[offset as usize] {
            self.edges[edge].multiplicity += 1;
            self.edges[edge].current_sample += 1;
            self.increment_backwards(sequence, previous, offset - 1);
        }
    }

    fn thread(&mut self, sequence: &[u8], start: usize, stop: usize, is_ref: bool) {
        let start_pos = if is_ref {
            Some(start)
        } else {
            (start..stop.saturating_sub(self.kmer)).find(|position| {
                !self
                    .non_unique
                    .contains(&sequence[*position..*position + self.kmer])
            })
        };
        let Some(start_pos) = start_pos else {
            return;
        };
        let first_kmer = &sequence[start_pos..start_pos + self.kmer];
        let mut vertex = self.vertex_for(first_kmer, true);
        self.increment_backwards(first_kmer, vertex, self.kmer as isize - 2);
        if is_ref {
            self.ref_source = vertex;
            self.ref_source_kmer = sequence[start..start + self.kmer].to_vec();
        }
        for position in start_pos + 1..=stop - self.kmer {
            let suffix = sequence[position + self.kmer - 1];
            let existing = self.live_outgoing_iter(vertex).find(|edge| {
                let target = self.edges[*edge].target;
                self.vertices[target].last == suffix
            });
            if let Some(edge) = existing {
                self.edges[edge].multiplicity += 1;
                self.edges[edge].current_sample += 1;
                vertex = self.edges[edge].target;
            } else {
                let kmer = &sequence[position..position + self.kmer];
                let next = self.vertex_for(kmer, false);
                self.add_edge(vertex, next, is_ref, 1);
                vertex = next;
            }
        }
        if is_ref {
            self.ref_sink = vertex;
        }
    }

    fn remove_edge(&mut self, edge: EdgeId) {
        if !self.edges[edge].alive {
            return;
        }
        self.edges[edge].alive = false;
        let source = self.edges[edge].source;
        let target = self.edges[edge].target;
        self.vertices[source].live_outgoing -= 1;
        self.vertices[target].live_incoming -= 1;
    }

    fn remove_vertex(&mut self, vertex: VertexId) {
        if !self.vertices[vertex].alive {
            return;
        }
        self.vertices[vertex].alive = false;
        self.live_vertices -= 1;
        // ReadThreadingGraph.removeVertex also drops the vertex's kmer from
        // kmerToVertexMap, which isLowQualityGraph counts after pruning.
        if self.track_unique {
            let sequence = std::mem::take(&mut self.vertices[vertex].sequence);
            self.unique.remove(&sequence);
            self.vertices[vertex].sequence = sequence;
        }
        // remove_edge never touches the adjacency lists themselves.
        for index in 0..self.vertices[vertex].incoming.len() {
            let edge = self.vertices[vertex].incoming[index];
            self.remove_edge(edge);
        }
        for index in 0..self.vertices[vertex].outgoing.len() {
            let edge = self.vertices[vertex].outgoing[index];
            self.remove_edge(edge);
        }
    }

    fn sources(&self) -> Vec<VertexId> {
        self.vertices
            .iter()
            .enumerate()
            .filter(|(_, vertex)| vertex.alive && vertex.live_incoming == 0)
            .map(|(id, _)| id)
            .collect()
    }

    fn prune(&mut self) {
        let sources = self.sources();
        let mut starts: VecDeque<VertexId> = sources.iter().copied().collect();
        let mut seen = vec![false; self.vertices.len()];
        for source in sources {
            seen[source] = true;
        }
        let mut chains: Vec<Vec<EdgeId>> = Vec::new();
        while let Some(start) = starts.pop_front() {
            for first in self.live_outgoing(start) {
                let mut chain = vec![first];
                let mut last = self.edges[first].target;
                loop {
                    let outgoing = self.live_outgoing(last);
                    if self.vertices[last].live_outgoing != 1
                        || self.vertices[last].live_incoming > 1
                        || last == start
                    {
                        break;
                    }
                    chain.push(outgoing[0]);
                    last = self.edges[outgoing[0]].target;
                }
                chains.push(chain);
                if !seen[last] {
                    seen[last] = true;
                    starts.push_back(last);
                }
            }
        }
        for chain in chains {
            if chain
                .iter()
                .all(|edge| self.edges[*edge].pruning_multiplicity < 2 && !self.edges[*edge].is_ref)
            {
                for edge in chain {
                    self.remove_edge(edge);
                }
            }
        }
        let orphans: Vec<_> = self
            .vertices
            .iter()
            .enumerate()
            .filter(|(id, vertex)| {
                vertex.alive
                    && *id != self.ref_source
                    && vertex.live_incoming == 0
                    && vertex.live_outgoing == 0
            })
            .map(|(id, _)| id)
            .collect();
        for vertex in orphans {
            self.remove_vertex(vertex);
        }
    }

    fn has_cycle(&self) -> bool {
        fn visit(graph: &Graph, vertex: VertexId, color: &mut [u8]) -> bool {
            color[vertex] = 1;
            for edge in graph.live_outgoing_iter(vertex) {
                let target = graph.edges[edge].target;
                if color[target] == 1 || (color[target] == 0 && visit(graph, target, color)) {
                    return true;
                }
            }
            color[vertex] = 2;
            false
        }
        let mut color = vec![0_u8; self.vertices.len()];
        self.vertices
            .iter()
            .enumerate()
            .any(|(id, vertex)| vertex.alive && color[id] == 0 && visit(self, id, &mut color))
    }

    fn reachable(&self, origin: VertexId, forward: bool) -> Vec<bool> {
        let mut found = vec![false; self.vertices.len()];
        let mut queue = VecDeque::from([origin]);
        while let Some(vertex) = queue.pop_front() {
            if found[vertex] {
                continue;
            }
            found[vertex] = true;
            if forward {
                for edge in self.live_outgoing_iter(vertex) {
                    queue.push_back(self.edges[edge].target);
                }
            } else {
                for edge in self.live_incoming_iter(vertex) {
                    queue.push_back(self.edges[edge].source);
                }
            }
        }
        found
    }

    fn retain_reference_paths(&mut self) {
        let from_source = self.reachable(self.ref_source, true);
        let to_sink = self.reachable(self.ref_sink, false);
        let remove: Vec<_> = self
            .vertices
            .iter()
            .enumerate()
            .filter(|(id, vertex)| vertex.alive && !(from_source[*id] && to_sink[*id]))
            .map(|(id, _)| id)
            .collect();
        for vertex in remove {
            self.remove_vertex(vertex);
        }
    }

    fn low_quality(&self) -> bool {
        self.non_unique.len() * 4 > self.unique.len()
    }

    fn alive_vertices(&self) -> Vec<VertexId> {
        self.vertices
            .iter()
            .enumerate()
            .filter(|(_, vertex)| vertex.alive)
            .map(|(id, _)| id)
            .collect()
    }

    fn is_reference_node(&self, vertex: VertexId) -> bool {
        self.live_incoming_iter(vertex)
            .chain(self.live_outgoing_iter(vertex))
            .any(|edge| self.edges[edge].is_ref)
            || self.live_vertices == 1
    }

    fn reference_source(&self) -> Option<VertexId> {
        self.alive_vertices().into_iter().find(|vertex| {
            !self
                .live_incoming_iter(*vertex)
                .any(|edge| self.edges[edge].is_ref)
                && self
                    .live_outgoing_iter(*vertex)
                    .any(|edge| self.edges[edge].is_ref)
        })
    }

    fn reference_sink(&self) -> Option<VertexId> {
        self.alive_vertices().into_iter().find(|vertex| {
            !self
                .live_outgoing_iter(*vertex)
                .any(|edge| self.edges[edge].is_ref)
                && self
                    .live_incoming_iter(*vertex)
                    .any(|edge| self.edges[edge].is_ref)
        })
    }

    fn outgoing_vertices(&self, vertex: VertexId) -> Vec<VertexId> {
        let mut result = Vec::new();
        for edge in self.live_outgoing_iter(vertex) {
            let target = self.edges[edge].target;
            if !result.contains(&target) {
                result.push(target);
            }
        }
        result
    }

    fn incoming_vertices(&self, vertex: VertexId) -> Vec<VertexId> {
        let mut result = Vec::new();
        for edge in self.live_incoming_iter(vertex) {
            let source = self.edges[edge].source;
            if !result.contains(&source) {
                result.push(source);
            }
        }
        result
    }

    fn convert_to_sequence_graph(&mut self) {
        // Only sequences change below, so "source" is the live in-degree now.
        for vertex in 0..self.vertices.len() {
            if self.vertices[vertex].alive && self.vertices[vertex].live_incoming != 0 {
                let suffix = *self.vertices[vertex]
                    .sequence
                    .last()
                    .expect("read-threading vertex is nonempty");
                self.vertices[vertex].sequence = vec![suffix];
            }
        }
        for edge in &mut self.edges {
            edge.current_sample = 0;
            edge.pruning_multiplicity = edge.multiplicity;
        }
    }

    fn zip_linear_chains(&mut self) -> bool {
        let starts: Vec<_> = self
            .alive_vertices()
            .into_iter()
            .filter(|vertex| {
                if self.vertices[*vertex].live_outgoing != 1 {
                    return false;
                }
                if self.vertices[*vertex].live_incoming != 1 {
                    return true;
                }
                let incoming = self
                    .live_incoming_iter(*vertex)
                    .next()
                    .expect("live incoming edge");
                self.vertices[self.edges[incoming].source].live_outgoing > 1
            })
            .collect();
        if starts.is_empty() {
            return false;
        }
        let mut merged = false;
        for start in starts {
            if !self.vertices[start].alive {
                continue;
            }
            let mut chain = vec![start];
            let mut last = start;
            let mut last_ref = self.is_reference_node(last);
            loop {
                if self.vertices[last].live_outgoing != 1 {
                    break;
                }
                let outgoing = self
                    .live_outgoing_iter(last)
                    .next()
                    .expect("live outgoing edge");
                let target = self.edges[outgoing].target;
                if self.vertices[target].live_incoming != 1 || target == last {
                    break;
                }
                let target_ref = self.is_reference_node(target);
                if last_ref != target_ref {
                    break;
                }
                chain.push(target);
                last = target;
                last_ref = target_ref;
            }
            if chain.len() == 1 {
                continue;
            }
            let mut sequence = Vec::new();
            for vertex in &chain {
                sequence.extend_from_slice(&self.vertices[*vertex].sequence);
            }
            let new_vertex = self.add_vertex(sequence);
            for edge in self.live_outgoing(last) {
                let old = self.edges[edge].clone();
                self.add_base_edge(new_vertex, old.target, old.is_ref, old.multiplicity);
            }
            for edge in self.live_incoming(start) {
                let old = self.edges[edge].clone();
                self.add_base_edge(old.source, new_vertex, old.is_ref, old.multiplicity);
            }
            for vertex in chain {
                self.remove_vertex(vertex);
            }
            merged = true;
        }
        merged
    }

    fn common_prefix_suffix(&self, vertices: &[VertexId]) -> (Vec<u8>, Vec<u8>) {
        let first = &self.vertices[vertices[0]].sequence;
        let min_len = vertices
            .iter()
            .map(|vertex| self.vertices[*vertex].sequence.len())
            .min()
            .unwrap_or(0);
        let mut prefix = 0;
        while prefix < min_len
            && vertices
                .iter()
                .all(|vertex| self.vertices[*vertex].sequence[prefix] == first[prefix])
        {
            prefix += 1;
        }
        let mut suffix = 0;
        while suffix < min_len - prefix
            && vertices.iter().all(|vertex| {
                let sequence = &self.vertices[*vertex].sequence;
                sequence[sequence.len() - 1 - suffix] == first[first.len() - 1 - suffix]
            })
        {
            suffix += 1;
        }
        (
            first[..prefix].to_vec(),
            first[first.len() - suffix..].to_vec(),
        )
    }

    fn split_shared_vertices(
        &mut self,
        middles: &[VertexId],
        top: Option<VertexId>,
        bottom: Option<VertexId>,
    ) {
        #[derive(Clone)]
        struct TempEdge {
            source: usize,
            target: usize,
            is_ref: bool,
            multiplicity: i32,
        }
        const PREFIX: usize = usize::MAX - 1;
        const SUFFIX: usize = usize::MAX;
        let (prefix_sequence, suffix_sequence) = self.common_prefix_suffix(middles);
        let mut remaining = Vec::new();
        let mut temp_edges: Vec<TempEdge> = Vec::new();
        for middle in middles {
            let incoming = self.live_incoming(*middle);
            let outgoing = self.live_outgoing(*middle);
            let to_mid = incoming
                .first()
                .map(|edge| self.edges[*edge].clone())
                .unwrap_or(Edge {
                    source: 0,
                    target: 0,
                    multiplicity: 0,
                    current_sample: 0,
                    pruning_multiplicity: 0,
                    is_ref: self.is_reference_node(*middle),
                    alive: true,
                });
            let from_mid = outgoing
                .first()
                .map(|edge| self.edges[*edge].clone())
                .unwrap_or(Edge {
                    source: 0,
                    target: 0,
                    multiplicity: 0,
                    current_sample: 0,
                    pruning_multiplicity: 0,
                    is_ref: self.is_reference_node(*middle),
                    alive: true,
                });
            let sequence = &self.vertices[*middle].sequence;
            let middle_end = sequence.len() - suffix_sequence.len();
            let center = &sequence[prefix_sequence.len()..middle_end];
            if !center.is_empty() {
                let temp_id = remaining.len();
                remaining.push(center.to_vec());
                temp_edges.push(TempEdge {
                    source: PREFIX,
                    target: temp_id,
                    is_ref: to_mid.is_ref,
                    multiplicity: to_mid.multiplicity,
                });
                temp_edges.push(TempEdge {
                    source: temp_id,
                    target: SUFFIX,
                    is_ref: from_mid.is_ref,
                    multiplicity: from_mid.multiplicity,
                });
            } else if let Some(edge) = temp_edges
                .iter_mut()
                .find(|edge| edge.source == PREFIX && edge.target == SUFFIX)
            {
                edge.multiplicity += to_mid.multiplicity + from_mid.multiplicity;
                edge.is_ref |= to_mid.is_ref || from_mid.is_ref;
            } else {
                temp_edges.push(TempEdge {
                    source: PREFIX,
                    target: SUFFIX,
                    is_ref: to_mid.is_ref || from_mid.is_ref,
                    multiplicity: to_mid.multiplicity + from_mid.multiplicity,
                });
            }
        }
        for middle in middles {
            self.remove_vertex(*middle);
        }
        let remaining_ids: Vec<_> = remaining
            .into_iter()
            .map(|sequence| self.add_vertex(sequence))
            .collect();
        let prefix_suffix_edge = temp_edges
            .iter()
            .any(|edge| edge.source == PREFIX && edge.target == SUFFIX);
        let prefix_out_degree = temp_edges
            .iter()
            .filter(|edge| edge.source == PREFIX)
            .count();
        let only_prefix_suffix = prefix_suffix_edge && prefix_out_degree == 1;
        let need_prefix = !prefix_sequence.is_empty() || (top.is_none() && !only_prefix_suffix);
        let need_suffix = !suffix_sequence.is_empty() || (bottom.is_none() && !only_prefix_suffix);
        let prefix_id = need_prefix.then(|| self.add_vertex(prefix_sequence));
        let suffix_id = need_suffix.then(|| self.add_vertex(suffix_sequence));
        let top_connect = prefix_id.or(top);
        let bottom_connect = suffix_id.or(bottom);
        if let (Some(top_vertex), Some(prefix_vertex)) = (top, prefix_id) {
            let is_ref = temp_edges
                .iter()
                .filter(|edge| edge.source == PREFIX)
                .any(|edge| edge.is_ref);
            self.add_base_edge(top_vertex, prefix_vertex, is_ref, 1);
        }
        if let (Some(suffix_vertex), Some(bottom_vertex)) = (suffix_id, bottom) {
            let is_ref = temp_edges
                .iter()
                .filter(|edge| edge.target == SUFFIX)
                .any(|edge| edge.is_ref);
            self.add_base_edge(suffix_vertex, bottom_vertex, is_ref, 1);
        }
        if let Some(top_vertex) = top_connect {
            for edge in temp_edges.iter().filter(|edge| edge.source == PREFIX) {
                let target = if edge.target == SUFFIX {
                    bottom_connect
                } else {
                    Some(remaining_ids[edge.target])
                };
                if let Some(target) = target {
                    self.add_base_edge(top_vertex, target, edge.is_ref, edge.multiplicity);
                }
            }
        }
        if let Some(bottom_vertex) = bottom_connect {
            for edge in temp_edges.iter().filter(|edge| edge.target == SUFFIX) {
                let source = if edge.source == PREFIX {
                    top_connect
                } else {
                    Some(remaining_ids[edge.source])
                };
                if edge.source != PREFIX
                    && let Some(source) = source
                {
                    self.add_base_edge(source, bottom_vertex, edge.is_ref, edge.multiplicity);
                }
            }
        }
    }

    fn transform_diamond(&mut self, top: VertexId) -> bool {
        if self.vertices[top].live_outgoing <= 1 {
            return false;
        }
        let middles = self.outgoing_vertices(top);
        if middles.len() <= 1 {
            return false;
        }
        let mut bottom = None;
        for middle in &middles {
            if self.vertices[*middle].live_outgoing == 0
                || self.vertices[*middle].live_incoming != 1
            {
                return false;
            }
            for target in self.outgoing_vertices(*middle) {
                if bottom.is_some_and(|known| known != target) {
                    return false;
                }
                bottom = Some(target);
            }
        }
        let Some(bottom) = bottom else { return false };
        if self.vertices[bottom].live_incoming != middles.len() {
            return false;
        }
        let (prefix, suffix) = self.common_prefix_suffix(&middles);
        if prefix.is_empty() && suffix.is_empty() {
            return false;
        }
        self.split_shared_vertices(&middles, Some(top), Some(bottom));
        true
    }

    fn transform_tail(&mut self, top: VertexId) -> bool {
        if self.vertices[top].live_outgoing <= 1 {
            return false;
        }
        let tails = self.outgoing_vertices(top);
        if tails.len() <= 1
            || tails.iter().any(|tail| {
                self.vertices[*tail].live_outgoing != 0 || self.vertices[*tail].live_incoming > 1
            })
        {
            return false;
        }
        let (_, suffix) = self.common_prefix_suffix(&tails);
        if suffix.len() < 10 {
            return false;
        }
        self.split_shared_vertices(&tails, Some(top), None);
        true
    }

    fn common_suffix_split(&mut self, bottom: VertexId) -> bool {
        if self.vertices[bottom].live_incoming < 2 {
            return false;
        }
        let middles = self.incoming_vertices(bottom);
        if middles.len() < 2 {
            return false;
        }
        let bottom_outgoing = self.outgoing_vertices(bottom);
        for middle in &middles {
            let outgoing = self.live_outgoing(*middle);
            if *middle == bottom
                || outgoing.len() != 1
                || self.edges[outgoing[0]].target != bottom
                || bottom_outgoing.contains(middle)
            {
                return false;
            }
        }
        let min_len = middles
            .iter()
            .map(|middle| self.vertices[*middle].sequence.len())
            .min()
            .unwrap_or(0);
        let first = &self.vertices[middles[0]].sequence;
        let mut suffix_len = 0;
        while suffix_len < min_len
            && middles.iter().all(|middle| {
                let sequence = &self.vertices[*middle].sequence;
                sequence[sequence.len() - 1 - suffix_len] == first[first.len() - 1 - suffix_len]
            })
        {
            suffix_len += 1;
        }
        if suffix_len == 0
            || middles
                .iter()
                .all(|middle| self.vertices[*middle].sequence.len() == suffix_len)
        {
            return false;
        }
        if let Some(ref_source) = self.reference_source()
            && middles.contains(&ref_source)
            && self.vertices[ref_source].sequence.len() == suffix_len
        {
            return false;
        }
        let suffix = first[first.len() - suffix_len..].to_vec();
        for middle in middles.clone() {
            let suffix_vertex = self.add_vertex(suffix.clone());
            let outgoing = self.live_outgoing(middle)[0];
            let out = self.edges[outgoing].clone();
            let prefix_len = self.vertices[middle].sequence.len() - suffix_len;
            let incoming_target = if prefix_len == 0 {
                suffix_vertex
            } else {
                let prefix = self.add_vertex(self.vertices[middle].sequence[..prefix_len].to_vec());
                self.add_base_edge(prefix, suffix_vertex, out.is_ref, 1);
                prefix
            };
            self.add_base_edge(suffix_vertex, out.target, out.is_ref, out.multiplicity);
            for edge in self.live_incoming(middle) {
                let old = self.edges[edge].clone();
                self.add_base_edge(old.source, incoming_target, old.is_ref, old.multiplicity);
            }
        }
        for middle in middles {
            self.remove_vertex(middle);
        }
        true
    }

    fn merge_common_suffix(&mut self, vertex: VertexId) -> bool {
        if self.vertices[vertex].live_incoming == 0 {
            return false;
        }
        let previous = self.incoming_vertices(vertex);
        if previous.is_empty() {
            return false;
        }
        let sequence = self.vertices[previous[0]].sequence.clone();
        if previous.iter().any(|prev| {
            self.vertices[*prev].sequence != sequence
                || self.outgoing_vertices(*prev) != vec![vertex]
                || self.vertices[*prev].live_incoming == 0
        }) {
            return false;
        }
        let mut merged_sequence = sequence;
        merged_sequence.extend_from_slice(&self.vertices[vertex].sequence);
        let merged = self.add_vertex(merged_sequence);
        for prev in &previous {
            for edge in self.live_incoming(*prev) {
                let old = self.edges[edge].clone();
                self.add_base_edge(old.source, merged, old.is_ref, old.multiplicity);
            }
        }
        for edge in self.live_outgoing(vertex) {
            let old = self.edges[edge].clone();
            self.add_base_edge(merged, old.target, old.is_ref, old.multiplicity);
        }
        for prev in previous {
            self.remove_vertex(prev);
        }
        self.remove_vertex(vertex);
        true
    }

    fn transform_until(&mut self, kind: u8) -> bool {
        let mut any = false;
        let mut already_split: Vec<bool> = Vec::new();
        loop {
            let mut changed = false;
            // A transform that returns false changes nothing, and the scan
            // stops at the first one that returns true, so walking the live
            // vertex list equals walking a snapshot of it.
            for vertex in 0..self.vertices.len() {
                if !self.vertices[vertex].alive {
                    continue;
                }
                let transformed = match kind {
                    0 => self.transform_diamond(vertex),
                    1 => self.transform_tail(vertex),
                    2 => {
                        if already_split.len() <= vertex {
                            already_split.resize(self.vertices.len().max(vertex + 1), false);
                        }
                        if !already_split[vertex] {
                            already_split[vertex] = true;
                            self.common_suffix_split(vertex)
                        } else {
                            false
                        }
                    }
                    3 => self.merge_common_suffix(vertex),
                    _ => unreachable!(),
                };
                if transformed {
                    any = true;
                    changed = true;
                    break;
                }
            }
            if !changed {
                break;
            }
        }
        any
    }

    fn simplify(&mut self) {
        self.zip_linear_chains();
        for _ in 0..100 {
            let mut changed = false;
            changed |= self.transform_until(0);
            changed |= self.transform_until(1);
            changed |= self.transform_until(2);
            changed |= self.transform_until(3);
            changed |= self.zip_linear_chains();
            if !changed {
                break;
            }
        }
    }

    fn cleanup_sequence_graph(&mut self) {
        self.zip_linear_chains();
        self.simplify();
        if let (Some(source), Some(sink)) = (self.reference_source(), self.reference_sink()) {
            self.ref_source = source;
            self.ref_sink = sink;
            self.retain_reference_paths();
            self.simplify();
            if let Some(source) = self.reference_source() {
                self.ref_source = source;
            }
            if let Some(sink) = self.reference_sink() {
                self.ref_sink = sink;
            }
        }
    }

    fn heaviest_incoming(&self, vertex: VertexId) -> Option<EdgeId> {
        let mut edges = self.live_incoming_iter(vertex);
        let mut best = edges.next()?;
        for edge in edges {
            if self.edges[edge].multiplicity > self.edges[best].multiplicity {
                best = edge;
            }
        }
        Some(best)
    }

    fn next_reference_vertex(
        &self,
        vertex: VertexId,
        blacklist: Option<EdgeId>,
    ) -> Option<VertexId> {
        if let Some(edge) = self
            .live_outgoing_iter(vertex)
            .find(|edge| self.edges[*edge].is_ref)
        {
            return Some(self.edges[edge].target);
        }
        let allowed: Vec<_> = self
            .live_outgoing_iter(vertex)
            .filter(|edge| Some(*edge) != blacklist)
            .take(2)
            .collect();
        (allowed.len() == 1).then(|| self.edges[allowed[0]].target)
    }

    fn previous_reference_vertex(&self, vertex: VertexId) -> Option<VertexId> {
        self.live_incoming_iter(vertex)
            .map(|edge| self.edges[edge].source)
            .find(|source| self.is_reference_node(*source))
    }

    fn reference_path(
        &self,
        start: VertexId,
        down: bool,
        blacklist: Option<EdgeId>,
    ) -> Vec<VertexId> {
        let mut result = Vec::new();
        let mut current = Some(start);
        while let Some(vertex) = current {
            result.push(vertex);
            current = if down {
                self.next_reference_vertex(vertex, blacklist)
            } else {
                self.previous_reference_vertex(vertex)
            };
        }
        result
    }

    fn bases_for_path(&self, path: &[VertexId], expand_source: bool) -> Vec<u8> {
        let mut result = Vec::new();
        for vertex in path {
            if expand_source && self.vertices[*vertex].live_incoming == 0 {
                result.extend(self.vertices[*vertex].sequence.iter().rev());
            } else {
                result.push(
                    *self.vertices[*vertex]
                        .sequence
                        .last()
                        .expect("nonempty kmer"),
                );
            }
        }
        result
    }

    fn find_tail_path(&self, start: VertexId) -> Option<Vec<VertexId>> {
        let mut path = VecDeque::new();
        let mut visited = FastSet::default();
        let mut vertex = start;
        while self.vertices[vertex].live_incoming == 1 && self.vertices[vertex].live_outgoing < 2 {
            let edge = self
                .live_incoming_iter(vertex)
                .next()
                .expect("live incoming edge");
            if self.edges[edge].pruning_multiplicity < 2 {
                visited.extend(path.drain(..));
            } else {
                path.push_front(vertex);
            }
            vertex = self.edges[edge].source;
            if path.contains(&vertex) || visited.contains(&vertex) {
                return None;
            }
        }
        path.push_front(vertex);
        (self.vertices[vertex].live_outgoing > 1).then(|| path.into_iter().collect())
    }

    fn find_head_path(&self, start: VertexId) -> Option<Vec<VertexId>> {
        let mut path = VecDeque::new();
        let mut visited = FastSet::default();
        let mut vertex = start;
        while !self.is_reference_node(vertex) && self.vertices[vertex].live_outgoing == 1 {
            let edge = self
                .live_outgoing_iter(vertex)
                .next()
                .expect("live outgoing edge");
            if self.edges[edge].pruning_multiplicity < 2 {
                visited.extend(path.drain(..));
            } else {
                path.push_front(vertex);
            }
            vertex = self.edges[edge].target;
            if path.contains(&vertex) || visited.contains(&vertex) {
                return None;
            }
        }
        path.push_front(vertex);
        self.is_reference_node(vertex)
            .then(|| path.into_iter().collect())
    }

    fn cigar_elements(cigar: &[u32]) -> Vec<(usize, u8)> {
        cigar
            .iter()
            .map(|value| ((value >> 4) as usize, (value & 15) as u8))
            .collect()
    }

    fn recover_dangling_tail(&mut self, sink: VertexId, ref_source: Option<VertexId>) {
        let Some(alt_path) = self.find_tail_path(sink) else {
            return;
        };
        if alt_path.len() < 5 || ref_source == Some(alt_path[0]) {
            return;
        }
        let blacklist = self.heaviest_incoming(alt_path[1]);
        let ref_path = self.reference_path(alt_path[0], true, blacklist);
        let ref_bases = self.bases_for_path(&ref_path, false);
        let alt_bases = self.bases_for_path(&alt_path, false);
        let cigar = dangling_cigar(&ref_bases, &alt_bases);
        let elements = Self::cigar_elements(&cigar);
        if elements.is_empty()
            || elements.len() > 3
            || elements.last().is_none_or(|elem| elem.1 != 0)
        {
            return;
        }
        let reference_len: usize = elements
            .iter()
            .filter(|elem| matches!(elem.1, 0 | 2 | 3 | 7 | 8))
            .map(|elem| elem.0)
            .sum();
        let read_len: usize = elements
            .iter()
            .filter(|elem| matches!(elem.1, 0 | 1 | 4 | 7 | 8))
            .map(|elem| elem.0)
            .sum();
        let last_ref = reference_len - 1;
        let mut matching = 0;
        for length in 1..=alt_bases.len() {
            let Some(ref_index) = last_ref.checked_sub(length - 1) else {
                break;
            };
            if ref_bases[ref_index] != alt_bases[alt_bases.len() - length] {
                break;
            }
            matching = length;
        }
        matching = matching.min(elements.last().expect("nonempty").0);
        if matching == 0 {
            return;
        }
        let alt_index = read_len.saturating_sub(matching + 1);
        let leading_deletion = elements[0].1 == 2 && elements[0].0 + matching == last_ref + 1;
        let ref_index = last_ref - matching + 1 + usize::from(leading_deletion);
        if ref_index == 0 || alt_index >= alt_path.len() || ref_index >= ref_path.len() {
            return;
        }
        self.add_edge(alt_path[alt_index], ref_path[ref_index], false, 1);
    }

    fn extend_dangling_head(
        &mut self,
        alt_path: &mut Vec<VertexId>,
        ref_path: &[VertexId],
        count: usize,
        cigar: &[(usize, u8)],
    ) -> bool {
        let last_index = alt_path.len() - 1;
        let offset: isize = cigar
            .iter()
            .map(|(length, op)| {
                let ref_len = if matches!(op, 0 | 2 | 3 | 7 | 8) {
                    *length as isize
                } else {
                    0
                };
                let read_len = if matches!(op, 0 | 1 | 4 | 7 | 8) {
                    *length as isize
                } else {
                    0
                };
                ref_len - read_len
            })
            .sum();
        let ref_index = last_index as isize + offset + count as isize;
        if ref_index < 0 || ref_index as usize >= ref_path.len() {
            return false;
        }
        let dangling_source = alt_path.pop().expect("head path nonempty");
        let mut extended = self.vertices[ref_path[ref_index as usize]].sequence[..count].to_vec();
        extended.extend_from_slice(&self.vertices[dangling_source].sequence);
        let mut source_edges = self.live_outgoing_iter(dangling_source);
        let Some(mut source_edge) = source_edges.next() else {
            return false;
        };
        for edge in source_edges {
            if self.edges[edge].multiplicity > self.edges[source_edge].multiplicity {
                source_edge = edge;
            }
        }
        let multiplicity = self.edges[source_edge].multiplicity;
        let mut previous = self.edges[source_edge].target;
        self.remove_edge(source_edge);
        for index in (1..=count).rev() {
            let vertex = self.add_vertex(extended[index..index + self.kmer].to_vec());
            let edge = self.add_edge(vertex, previous, false, 1);
            self.edges[edge].multiplicity = multiplicity;
            alt_path.push(vertex);
            previous = vertex;
        }
        true
    }

    fn recover_dangling_head(&mut self, source: VertexId, ref_sink: Option<VertexId>) {
        let Some(mut alt_path) = self.find_head_path(source) else {
            return;
        };
        if alt_path.len() < 5 || ref_sink == Some(alt_path[0]) {
            return;
        }
        let ref_path = self.reference_path(alt_path[0], false, None);
        let ref_bases = self.bases_for_path(&ref_path, true);
        let alt_bases = self.bases_for_path(&alt_path, true);
        let cigar = dangling_cigar(&ref_bases, &alt_bases);
        let elements = Self::cigar_elements(&cigar);
        if elements.is_empty() || elements.len() > 3 || elements[0].1 != 0 {
            return;
        }
        let max_index = elements[0].0;
        if max_index > ref_bases.len() || max_index > alt_bases.len() {
            return;
        }
        let max_mismatches = (max_index / self.kmer).max(1);
        let mut mismatches = 0;
        let mut last_mismatch = None;
        for index in 0..max_index {
            if ref_bases[index] != alt_bases[index] {
                mismatches += 1;
                if mismatches > max_mismatches {
                    return;
                }
                last_mismatch = Some(index);
            }
        }
        let Some(index) = last_mismatch else { return };
        if index == 0 || index >= ref_path.len() - 1 {
            return;
        }
        if index >= alt_path.len() {
            let count = index - alt_path.len() + 2;
            if !self.extend_dangling_head(&mut alt_path, &ref_path, count, &elements) {
                return;
            }
        }
        self.add_edge(ref_path[index + 1], alt_path[index], false, 1);
    }

    fn recover_dangling_branches(&mut self) {
        // Recovery adds only non-reference edges and vertices, and removes
        // only an edge leaving a non-reference vertex, so the reference
        // source and sink are fixed for the whole pass.
        let ref_source = self.reference_source();
        let ref_sink = self.reference_sink();
        let vertices = self.alive_vertices();
        for vertex in vertices {
            if self.vertices[vertex].alive
                && self.vertices[vertex].live_outgoing == 0
                && ref_sink != Some(vertex)
            {
                self.recover_dangling_tail(vertex, ref_source);
            }
        }
        let heads: Vec<_> = self
            .alive_vertices()
            .into_iter()
            .filter(|vertex| {
                self.vertices[*vertex].live_incoming == 0 && ref_source != Some(*vertex)
            })
            .collect();
        for vertex in heads {
            self.recover_dangling_head(vertex, ref_sink);
        }
    }
}

fn non_unique_kmers(sequences: &[(&[u8], usize, usize)], kmer: usize) -> KmerSet {
    let mut result = KmerSet::new(kmer);
    // One membership table, cleared per sequence: clear() keeps capacity.
    let longest = sequences
        .iter()
        .map(|(_, _, stop)| *stop)
        .max()
        .unwrap_or(0);
    let mut in_sequence = KmerSet::with_capacity(kmer, longest);
    for (sequence, _start, stop) in sequences {
        in_sequence.clear();
        // GATK's determineNonUniqueKmers intentionally (and surprisingly)
        // ignores SequenceForKmers.start here and scans from zero through
        // stop-k.  This matters for a read split after a low-quality base.
        if kmer <= 32 {
            let mask = if kmer == 32 {
                u64::MAX
            } else {
                (1_u64 << (kmer * 2)) - 1
            };
            let mut packed = 0_u64;
            let mut canonical = 0_usize;
            for (index, base) in sequence[..*stop].iter().copied().enumerate() {
                if let Some(bits) = base_bits(base) {
                    packed = ((packed << 2) | bits) & mask;
                    canonical += 1;
                } else {
                    packed = 0;
                    canonical = 0;
                }
                if index + 1 < kmer {
                    continue;
                }
                let start = index + 1 - kmer;
                let duplicate = if canonical >= kmer {
                    !in_sequence.insert_packed(packed)
                } else {
                    !in_sequence.insert(&sequence[start..=index])
                };
                if duplicate {
                    if canonical >= kmer {
                        result.insert_packed(packed);
                    } else {
                        result.insert(&sequence[start..=index]);
                    }
                }
            }
        } else {
            for position in 0..=stop.saturating_sub(kmer) {
                let current = &sequence[position..position + kmer];
                if !in_sequence.insert(current) {
                    result.insert(current);
                }
            }
        }
    }
    result
}

fn ref_has_repeated_kmer(reference: &[u8], kmer: usize) -> bool {
    if reference.len() < kmer {
        return false;
    }
    let mut seen = KmerSet::new(kmer);
    reference.windows(kmer).any(|current| !seen.insert(current))
}

#[derive(Debug)]
pub(super) struct PathResult {
    pub(super) bases: Vec<u8>,
    pub(super) score: f64,
}

#[derive(Clone)]
struct Candidate {
    vertex: VertexId,
    bases: Vec<u8>,
    score: f64,
}

impl Eq for Candidate {}

impl PartialEq for Candidate {
    fn eq(&self, other: &Self) -> bool {
        self.score.to_bits() == other.score.to_bits() && self.bases == other.bases
    }
}

impl Ord for Candidate {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score
            .total_cmp(&other.score)
            .then_with(|| self.bases.cmp(&other.bases))
    }
}

impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn k_best(graph: &Graph, limit: usize) -> Vec<PathResult> {
    let mut queue = BinaryHeap::new();
    queue.push(Candidate {
        vertex: graph.ref_source,
        bases: graph.vertices[graph.ref_source].sequence.clone(),
        score: 0.0,
    });
    let mut expanded = vec![0_usize; graph.vertices.len()];
    let mut result = Vec::new();
    while let Some(path) = queue.pop() {
        if path.vertex == graph.ref_sink {
            result.push(PathResult {
                bases: path.bases,
                score: path.score,
            });
            if result.len() == limit {
                break;
            }
            continue;
        }
        if expanded[path.vertex] >= limit {
            continue;
        }
        expanded[path.vertex] += 1;
        let total: i32 = graph
            .live_outgoing_iter(path.vertex)
            .map(|edge| graph.edges[edge].multiplicity)
            .sum();
        for edge in graph.live_outgoing_iter(path.vertex) {
            let target = graph.edges[edge].target;
            let mut bases = path.bases.clone();
            bases.extend_from_slice(&graph.vertices[target].sequence);
            let penalty = jmath::log10(f64::from(graph.edges[edge].multiplicity))
                - jmath::log10(f64::from(total));
            let score = path.score + penalty;
            queue.push(Candidate {
                vertex: target,
                bases,
                score,
            });
        }
    }
    result
}

pub(super) fn assemble_graph(
    reference: &[u8],
    reads: &[(Vec<u8>, Vec<u8>)],
    kmer: usize,
    allow_low_complexity: bool,
) -> Option<Vec<PathResult>> {
    if reference.len() < kmer || (!allow_low_complexity && ref_has_repeated_kmer(reference, kmer)) {
        return None;
    }
    let mut pending: Vec<(&[u8], usize, usize)> = vec![(reference, 0, reference.len())];
    for (bases, qualities) in reads {
        let mut last_good = None;
        for end in 0..=bases.len() {
            let usable = end < bases.len() && bases[end] != b'N' && qualities[end] >= 10;
            if usable {
                last_good.get_or_insert(end);
            } else if let Some(start) = last_good.take()
                && end - start >= kmer
            {
                pending.push((bases, start, end));
            }
        }
    }
    let non_unique = non_unique_kmers(&pending, kmer);
    // Capacity hint only: most read kmers land on existing vertices.
    let mut graph = Graph::new(kmer, non_unique, reference.len() * 2);
    graph.thread(reference, 0, reference.len(), true);
    graph.flush_sample();
    for (sequence, start, stop) in pending.into_iter().skip(1) {
        graph.thread(sequence, start, stop, false);
    }
    graph.flush_sample();
    if super::debug_assembly()
        && let Some(path) = std::env::var_os("NOMICS_DEBUG_GRAPH_DUMP")
    {
        let mut text = String::new();
        for v in &graph.vertices {
            text.push_str(&String::from_utf8_lossy(&v.sequence));
            text.push('\n');
        }
        for e in &graph.edges {
            text.push_str(&format!(
                "E {} {} {}/{} {}\n",
                String::from_utf8_lossy(&graph.vertices[e.source].sequence),
                String::from_utf8_lossy(&graph.vertices[e.target].sequence),
                e.multiplicity,
                e.pruning_multiplicity,
                e.is_ref
            ));
        }
        let tag = String::from_utf8_lossy(&reference[..12.min(reference.len())]).into_owned();
        let _ = std::fs::write(format!("{}.{tag}.k{kmer}", path.to_string_lossy()), text);
    }
    graph.prune();
    let cycle = graph.has_cycle();
    let low = !allow_low_complexity && graph.low_quality();
    if super::debug_assembly() {
        eprintln!(
            "    k={kmer}: non-unique {} unique {} cycle={cycle} low_quality={low}",
            graph.non_unique.len(),
            graph.unique.len()
        );
    }
    if cycle || low {
        return None;
    }
    graph.track_unique = false;
    graph.recover_dangling_branches();
    graph.retain_reference_paths();
    graph.convert_to_sequence_graph();
    graph.cleanup_sequence_graph();
    if super::debug_assembly()
        && let Some(path) = std::env::var_os("NOMICS_DEBUG_GRAPH_DUMP")
    {
        let mut text = String::new();
        for (i, e) in graph.edges.iter().enumerate() {
            if !e.alive {
                continue;
            }
            text.push_str(&format!(
                "E{i} {}:{} -> {}:{} m={} ref={}\n",
                e.source,
                String::from_utf8_lossy(&graph.vertices[e.source].sequence),
                e.target,
                String::from_utf8_lossy(&graph.vertices[e.target].sequence),
                e.multiplicity,
                e.is_ref
            ));
        }
        let tag = String::from_utf8_lossy(&reference[..12.min(reference.len())]).into_owned();
        let _ = std::fs::write(format!("{}.{tag}.k{kmer}.final", path.to_string_lossy()), text);
    }
    Some(k_best(&graph, 128))
}
