//! HashMap iteration order for distinct Comparable keys, OpenJDK 17.
//! Replaying insertion preserves resize splits and tree-bin root movement.
//! This matters for count ties and the order of floating-point additions.
use std::cmp::Ordering;

#[derive(Clone, Default)]
struct Node {
    parent: Option<usize>,
    left: Option<usize>,
    right: Option<usize>,
    red: bool,
}
#[derive(Default)]
struct Bucket {
    chain: Vec<usize>,
    tree: bool,
}
struct Map<'a, K> {
    entries: &'a [(u32, K)],
    buckets: Vec<Bucket>,
    nodes: Vec<Node>,
}
pub(super) fn order<K: Ord>(entries: &[(u32, K)]) -> Vec<usize> {
    let mut m = Map {
        entries,
        buckets: (0..16).map(|_| Bucket::default()).collect(),
        nodes: vec![Node::default(); entries.len()],
    };
    for (i, entry) in entries.iter().enumerate() {
        let index = entry.0 as usize & (m.buckets.len() - 1);
        if m.buckets[index].tree {
            let mut root = m.buckets[index].chain[0];
            let parent = m.insert_tree(root, i);
            let at = m.buckets[index]
                .chain
                .iter()
                .position(|&v| v == parent)
                .unwrap()
                + 1;
            m.buckets[index].chain.insert(at, i);
            root = m.balance(root, i);
            m.front(index, root);
        } else {
            m.buckets[index].chain.push(i);
            if m.buckets[index].chain.len() >= 9 {
                if m.buckets.len() < 64 {
                    m.resize();
                } else {
                    m.treeify(index);
                }
            }
        }
        if i + 1 > m.buckets.len() * 3 / 4 {
            m.resize();
        }
    }
    m.buckets.into_iter().flat_map(|b| b.chain).collect()
}
impl<K: Ord> Map<'_, K> {
    fn compare(&self, a: usize, b: usize) -> Ordering {
        (self.entries[a].0 as i32)
            .cmp(&(self.entries[b].0 as i32))
            .then_with(|| self.entries[a].1.cmp(&self.entries[b].1))
    }
    fn insert_tree(&mut self, root: usize, i: usize) -> usize {
        self.nodes[i] = Node {
            red: true,
            ..Node::default()
        };
        let mut p = root;
        loop {
            let left = self.compare(i, p) == Ordering::Less;
            if let Some(child) = if left {
                self.nodes[p].left
            } else {
                self.nodes[p].right
            } {
                p = child;
            } else {
                self.nodes[i].parent = Some(p);
                if left {
                    self.nodes[p].left = Some(i);
                } else {
                    self.nodes[p].right = Some(i);
                }
                return p;
            }
        }
    }
    fn front(&mut self, index: usize, root: usize) {
        let c = &mut self.buckets[index].chain;
        let at = c.iter().position(|&v| v == root).unwrap();
        c.remove(at);
        c.insert(0, root);
    }
    fn treeify(&mut self, index: usize) {
        let chain = self.buckets[index].chain.clone();
        let Some(&first) = chain.first() else {
            return;
        };
        self.nodes[first] = Node::default();
        let mut root = first;
        for &i in &chain[1..] {
            self.insert_tree(root, i);
            root = self.balance(root, i);
        }
        self.buckets[index].tree = true;
        self.front(index, root);
    }
    fn resize(&mut self) {
        let old_len = self.buckets.len();
        let old = std::mem::replace(
            &mut self.buckets,
            (0..old_len * 2).map(|_| Bucket::default()).collect(),
        );
        for (index, b) in old.into_iter().enumerate() {
            for &i in &b.chain {
                let slot = index
                    + if self.entries[i].0 as usize & old_len == 0 {
                        0
                    } else {
                        old_len
                    };
                self.buckets[slot].chain.push(i);
            }
            if b.tree {
                let both = !self.buckets[index].chain.is_empty()
                    && !self.buckets[index + old_len].chain.is_empty();
                for slot in [index, index + old_len] {
                    if self.buckets[slot].chain.len() > 6 {
                        self.buckets[slot].tree = true;
                        if both {
                            self.treeify(slot);
                        }
                    }
                }
            }
        }
    }
    fn red(&self, n: Option<usize>) -> bool {
        n.is_some_and(|i| self.nodes[i].red)
    }
    fn rotate(&mut self, mut root: usize, p: usize, left: bool) -> usize {
        let r = if left {
            self.nodes[p].right.unwrap()
        } else {
            self.nodes[p].left.unwrap()
        };
        let inner = if left {
            self.nodes[r].left
        } else {
            self.nodes[r].right
        };
        if left {
            self.nodes[p].right = inner;
        } else {
            self.nodes[p].left = inner;
        }
        if let Some(i) = inner {
            self.nodes[i].parent = Some(p);
        }
        let pp = self.nodes[p].parent;
        self.nodes[r].parent = pp;
        if let Some(pp) = pp {
            if self.nodes[pp].left == Some(p) {
                self.nodes[pp].left = Some(r);
            } else {
                self.nodes[pp].right = Some(r);
            }
        } else {
            root = r;
            self.nodes[r].red = false;
        }
        if left {
            self.nodes[r].left = Some(p);
        } else {
            self.nodes[r].right = Some(p);
        }
        self.nodes[p].parent = Some(r);
        root
    }
    fn balance(&mut self, mut root: usize, mut x: usize) -> usize {
        loop {
            let Some(p) = self.nodes[x].parent else {
                self.nodes[x].red = false;
                return x;
            };
            if !self.nodes[p].red {
                return root;
            }
            let Some(g) = self.nodes[p].parent else {
                return root;
            };
            let left = self.nodes[g].left == Some(p);
            let uncle = if left {
                self.nodes[g].right
            } else {
                self.nodes[g].left
            };
            if self.red(uncle) {
                self.nodes[uncle.unwrap()].red = false;
                self.nodes[p].red = false;
                self.nodes[g].red = true;
                x = g;
            } else {
                if (if left {
                    self.nodes[p].right
                } else {
                    self.nodes[p].left
                }) == Some(x)
                {
                    root = self.rotate(root, p, left);
                    x = p;
                }
                let p = self.nodes[x].parent.unwrap();
                let g = self.nodes[p].parent.unwrap();
                self.nodes[p].red = false;
                self.nodes[g].red = true;
                root = self.rotate(root, g, !left);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn tree_bins_reorder_count_ties() {
        // All sixteen Aa/BB concatenations have the same String.hashCode.
        // Default HashMap<String,Integer> on OpenJDK 17, after resize/treeify.
        let entries: Vec<_> = (0..16)
            .map(|n| {
                let s: String = (0..4)
                    .rev()
                    .map(|bit| if n & (1 << bit) == 0 { "Aa" } else { "BB" })
                    .collect();
                let h = s
                    .bytes()
                    .fold(0u32, |h, c| h.wrapping_mul(31).wrapping_add(c as u32));
                (h ^ (h >> 16), s)
            })
            .collect();
        assert_eq!(
            super::order(&entries),
            [3, 0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        );
    }
}
