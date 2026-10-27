//! Safe, index-arena port of bwa 0.7.19's `kbtree.h` insertion/search path.
//! Node indices replace C pointers; comparisons and duplicate routing follow
//! `KBTREE_INIT` exactly.

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct KeyLocation {
    node: usize,
    key: usize,
}

#[derive(Debug)]
struct Node<T> {
    is_internal: bool,
    keys: Vec<T>,
    children: Vec<usize>,
}
impl<T> Node<T> {
    fn leaf() -> Self {
        Self {
            is_internal: false,
            keys: Vec::new(),
            children: Vec::new(),
        }
    }
    fn internal() -> Self {
        Self {
            is_internal: true,
            keys: Vec::new(),
            children: Vec::new(),
        }
    }
}

/// B-tree with the same degree and duplicate-key behaviour as klib.
pub struct KBTree<T, C> {
    cmp: C,
    t: usize,
    capacity: usize,
    root: usize,
    nodes: Vec<Node<T>>,
    n_keys: usize,
}

impl<T: Copy, C: Fn(&T, &T) -> i32> KBTree<T, C> {
    /// Port of `kb_init_##name`; `None` is C's `NULL` result when `t < 2`.
    pub fn new(size: usize, elem_size: usize, cmp: C) -> Option<Self> {
        // C's formula is ((size - 4 - sizeof(void*)) / (sizeof(void*) +
        // sizeof(key_t)) + 1) >> 1; the required ABI fixes pointers at eight.
        if size < 12 || elem_size == 0 {
            return None;
        }
        let t = ((size - 4 - 8) / (8 + elem_size) + 1) >> 1;
        if t < 2 {
            return None;
        }
        Some(Self {
            cmp,
            t,
            capacity: 2 * t - 1,
            root: 0,
            nodes: vec![Node::leaf()],
            n_keys: 0,
        })
    }
    pub fn len(&self) -> usize {
        self.n_keys
    }
    pub fn is_empty(&self) -> bool {
        self.n_keys == 0
    }
    pub fn degree(&self) -> usize {
        self.t
    }

    /// Port of `__kb_getp_aux_##name`.
    fn getp_aux(&self, node: usize, key: &T) -> (isize, i32) {
        let x = &self.nodes[node];
        if x.keys.is_empty() {
            return (-1, 0);
        }
        let (mut begin, mut end) = (0, x.keys.len());
        while begin < end {
            let mid = (begin + end) >> 1;
            if (self.cmp)(&x.keys[mid], key) < 0 {
                begin = mid + 1;
            } else {
                end = mid;
            }
        }
        if begin == x.keys.len() {
            return ((x.keys.len() - 1) as isize, 1);
        }
        let r = (self.cmp)(key, &x.keys[begin]);
        if r < 0 {
            (begin as isize - 1, r)
        } else {
            (begin as isize, r)
        }
    }

    /// Port of `kb_getp_##name`, returning a reference into the stored key.
    pub fn getp(&self, key: &T) -> Option<&T> {
        let mut node = self.root;
        loop {
            let (i, r) = self.getp_aux(node, key);
            if i >= 0 && r == 0 {
                return Some(&self.nodes[node].keys[i as usize]);
            }
            if !self.nodes[node].is_internal {
                return None;
            }
            node = self.nodes[node].children[(i + 1) as usize];
        }
    }
    /// Mutable `kb_getp`, for BWA's mutation of a returned chain.
    pub fn getp_mut(&mut self, key: &T) -> Option<&mut T> {
        let p = self.get_location(key)?;
        Some(&mut self.nodes[p.node].keys[p.key])
    }
    fn get_location(&self, key: &T) -> Option<KeyLocation> {
        let mut node = self.root;
        loop {
            let (i, r) = self.getp_aux(node, key);
            if i >= 0 && r == 0 {
                return Some(KeyLocation {
                    node,
                    key: i as usize,
                });
            }
            if !self.nodes[node].is_internal {
                return None;
            }
            node = self.nodes[node].children[(i + 1) as usize];
        }
    }

    /// Port of `kb_intervalp_##name`.
    pub fn intervalp(&self, key: &T) -> (Option<&T>, Option<&T>) {
        let (lower, upper) = self.interval_locations(key);
        (
            lower.map(|p| &self.nodes[p.node].keys[p.key]),
            upper.map(|p| &self.nodes[p.node].keys[p.key]),
        )
    }
    /// Exact `kb_intervalp` locations. Handles permit safe mutation when C's
    /// lower and upper are the same pointer for an equal key.
    pub fn interval_locations(&self, key: &T) -> (Option<KeyLocation>, Option<KeyLocation>) {
        let mut node = self.root;
        let (mut lower, mut upper) = (None, None);
        loop {
            let (i, r) = self.getp_aux(node, key);
            if i >= 0 && r == 0 {
                let p = KeyLocation {
                    node,
                    key: i as usize,
                };
                return (Some(p), Some(p));
            }
            if i >= 0 {
                lower = Some(KeyLocation {
                    node,
                    key: i as usize,
                });
            }
            if i < self.nodes[node].keys.len() as isize - 1 {
                upper = Some(KeyLocation {
                    node,
                    key: (i + 1) as usize,
                });
            }
            if !self.nodes[node].is_internal {
                return (lower, upper);
            }
            node = self.nodes[node].children[(i + 1) as usize];
        }
    }
    pub fn get_location_mut(&mut self, p: KeyLocation) -> &mut T {
        &mut self.nodes[p.node].keys[p.key]
    }

    /// Mutable lower result of `kb_intervalp`. This is the safe equivalent of
    /// BWA mutating the `lower` pointer; `upper` is intentionally not returned
    /// simultaneously because equal keys make both C pointers alias.
    pub fn interval_lower_mut(&mut self, key: &T) -> Option<&mut T> {
        let (lower, _) = self.interval_locations(key);
        let p = lower?;
        Some(&mut self.nodes[p.node].keys[p.key])
    }

    /// Port of `__kb_split_##name`.
    fn split(&mut self, parent: usize, i: usize, child: usize) {
        let mut z = if self.nodes[child].is_internal {
            Node::internal()
        } else {
            Node::leaf()
        };
        let median = self.nodes[child].keys[self.t - 1];
        z.keys.extend_from_slice(&self.nodes[child].keys[self.t..]);
        if self.nodes[child].is_internal {
            z.children
                .extend_from_slice(&self.nodes[child].children[self.t..]);
            self.nodes[child].children.truncate(self.t);
        }
        self.nodes[child].keys.truncate(self.t - 1);
        let z_index = self.nodes.len();
        self.nodes.push(z);
        self.nodes[parent].children.insert(i + 1, z_index);
        self.nodes[parent].keys.insert(i, median);
    }
    /// Port of `__kb_putp_aux_##name`.
    fn putp_aux(&mut self, node: usize, key: T) {
        if !self.nodes[node].is_internal {
            let (i, _) = self.getp_aux(node, &key);
            self.nodes[node].keys.insert((i + 1) as usize, key);
            return;
        }
        let (at, _) = self.getp_aux(node, &key);
        let mut i = (at + 1) as usize;
        let child = self.nodes[node].children[i];
        if self.nodes[child].keys.len() == self.capacity {
            self.split(node, i, child);
            if (self.cmp)(&key, &self.nodes[node].keys[i]) > 0 {
                i += 1;
            }
        }
        self.putp_aux(self.nodes[node].children[i], key);
    }
    /// Port of `kb_putp_##name`; duplicate keys are intentionally retained.
    pub fn putp(&mut self, key: T) {
        self.n_keys += 1;
        let root = self.root;
        if self.nodes[root].keys.len() == self.capacity {
            let new_root = self.nodes.len();
            let mut s = Node::internal();
            s.children.push(root);
            self.nodes.push(s);
            self.root = new_root;
            self.split(new_root, 0, root);
        }
        self.putp_aux(self.root, key);
    }
    /// Port of `__kb_traverse`, yielding its in-order visit sequence.
    pub fn traverse<F: FnMut(&T)>(&self, mut visit: F) {
        // `position` alternates child and key visits: child 0, key 0,
        // child 1, key 1, ... child n. This is `__kb_traverse`'s stack.
        let mut stack = vec![(self.root, 0_usize)];
        while let Some((node, position)) = stack.pop() {
            let x = &self.nodes[node];
            if position > 2 * x.keys.len() {
                continue;
            }
            stack.push((node, position + 1));
            if position & 1 == 0 {
                if x.is_internal {
                    stack.push((x.children[position >> 1], 0));
                }
            } else {
                visit(&x.keys[position >> 1]);
            }
        }
    }
}
