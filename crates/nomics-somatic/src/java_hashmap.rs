//! String-key iteration order of OpenJDK 17 java/util/HashMap.java.
//! The groupingBy collector calls computeIfAbsent, including head insertion,
//! stable resize splitting, treeification and root-to-front tree-bin iteration.
use std::cmp::Ordering;
#[derive(Default)]
struct Node {
    key: String,
    hash: i32,
    values: Vec<usize>,
    parent: Option<usize>,
    child: [Option<usize>; 2],
    red: bool,
}
#[derive(Default, Clone)]
struct Bin {
    chain: Vec<usize>,
    tree: bool,
}
fn hash(s: &str) -> i32 {
    let h = s
        .encode_utf16()
        .fold(0u32, |h, c| h.wrapping_mul(31).wrapping_add(u32::from(c)));
    (h ^ (h >> 16)) as i32
}
fn compare(a: &Node, b: &Node) -> Ordering {
    a.hash
        .cmp(&b.hash)
        .then_with(|| a.key.encode_utf16().cmp(b.key.encode_utf16()))
}
fn rotate(nodes: &mut [Node], mut root: usize, p: usize, left: bool) -> usize {
    let dir = usize::from(left);
    let r = nodes[p].child[dir].unwrap();
    nodes[p].child[dir] = nodes[r].child[1 - dir];
    if let Some(rl) = nodes[r].child[1 - dir] {
        nodes[rl].parent = Some(p);
    }
    nodes[r].parent = nodes[p].parent;
    if let Some(pp) = nodes[p].parent {
        let side = usize::from(nodes[pp].child[0] != Some(p));
        nodes[pp].child[side] = Some(r);
    } else {
        root = r;
        nodes[r].red = false;
    }
    nodes[r].child[1 - dir] = Some(p);
    nodes[p].parent = Some(r);
    root
}
fn balance(nodes: &mut [Node], mut root: usize, mut x: usize) -> usize {
    nodes[x].red = true;
    loop {
        let Some(mut p) = nodes[x].parent else {
            nodes[x].red = false;
            return x;
        };
        if !nodes[p].red {
            return root;
        }
        let Some(mut gp) = nodes[p].parent else {
            return root;
        };
        let side = usize::from(nodes[gp].child[0] != Some(p));
        let uncle = nodes[gp].child[1 - side];
        if let Some(uncle) = uncle.filter(|&u| nodes[u].red) {
            nodes[uncle].red = false;
            nodes[p].red = false;
            nodes[gp].red = true;
            x = gp;
        } else {
            if nodes[p].child[1 - side] == Some(x) {
                x = p;
                root = rotate(nodes, root, x, side == 0);
                p = nodes[x].parent.unwrap();
                gp = nodes[p].parent.unwrap();
            }
            nodes[p].red = false;
            nodes[gp].red = true;
            root = rotate(nodes, root, gp, side != 0);
        }
    }
}
fn attach(nodes: &mut [Node], root: usize, x: usize) -> usize {
    let mut p = root;
    loop {
        let side = usize::from(compare(&nodes[x], &nodes[p]) == Ordering::Greater);
        if let Some(next) = nodes[p].child[side] {
            p = next;
        } else {
            nodes[p].child[side] = Some(x);
            nodes[x].parent = Some(p);
            return p;
        }
    }
}
fn move_front(bin: &mut Bin, root: usize) {
    let i = bin.chain.iter().position(|&v| v == root).unwrap();
    bin.chain.remove(i);
    bin.chain.insert(0, root);
}
fn treeify(nodes: &mut [Node], bin: &mut Bin) {
    let mut root = None;
    for &x in &bin.chain {
        nodes[x].parent = None;
        nodes[x].child = [None; 2];
        if let Some(r) = root {
            attach(nodes, r, x);
            root = Some(balance(nodes, r, x));
        } else {
            nodes[x].red = false;
            root = Some(x);
        }
    }
    bin.tree = true;
    move_front(bin, root.unwrap());
}
fn resize(nodes: &mut [Node], bins: &mut Vec<Bin>) {
    let n = bins.len();
    let mut next = vec![Bin::default(); n * 2];
    for (i, bin) in bins.iter().enumerate() {
        for &x in &bin.chain {
            let dest = if nodes[x].hash as usize & n == 0 {
                i
            } else {
                i + n
            };
            next[dest].chain.push(x);
        }
        if bin.tree {
            let both = !next[i].chain.is_empty() && !next[i + n].chain.is_empty();
            for index in [i, i + n] {
                if next[index].chain.len() > 6 {
                    if both {
                        treeify(nodes, &mut next[index]);
                    } else {
                        next[index].tree = true;
                    }
                }
            }
        }
    }
    *bins = next;
}
pub(crate) fn group_indices<'a>(keys: impl IntoIterator<Item = &'a str>) -> Vec<Vec<usize>> {
    let mut nodes: Vec<Node> = Vec::new();
    let mut bins = vec![Bin::default(); 16];
    for (i, key) in keys.into_iter().enumerate() {
        if nodes.len() > bins.len() * 3 / 4 {
            resize(&mut nodes, &mut bins);
        }
        let h = hash(key);
        let b = h as usize & (bins.len() - 1);
        if let Some(&x) = bins[b].chain.iter().find(|&&x| nodes[x].key == key) {
            nodes[x].values.push(i);
            continue;
        }
        let x = nodes.len();
        nodes.push(Node {
            key: key.into(),
            hash: h,
            values: vec![i],
            ..Node::default()
        });
        if bins[b].tree {
            let root = bins[b].chain[0];
            let parent = attach(&mut nodes, root, x);
            let at = bins[b].chain.iter().position(|&v| v == parent).unwrap() + 1;
            bins[b].chain.insert(at, x);
            let root = balance(&mut nodes, root, x);
            move_front(&mut bins[b], root);
        } else {
            let old_count = bins[b].chain.len();
            bins[b].chain.insert(0, x);
            if old_count >= 7 {
                if bins.len() < 64 {
                    resize(&mut nodes, &mut bins);
                } else {
                    treeify(&mut nodes, &mut bins[b]);
                }
            }
        }
    }
    bins.into_iter()
        .flat_map(|b| b.chain)
        .map(|x| std::mem::take(&mut nodes[x].values))
        .collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn java_tree_bins_and_resize() {
        // HashMapOracle.java: 32 String.hashCode collisions, then enough keys
        // to resize/split tree bins, then another colliding key and repeats.
        let mut keys = Vec::new();
        for i in 0..32 {
            keys.push(
                (0..5)
                    .map(|bit| if i & (1 << bit) == 0 { "Aa" } else { "BB" })
                    .collect::<String>(),
            );
        }
        for i in 0..100 {
            keys.push(format!("contig{i}"));
        }
        keys.extend([
            "AaAaAaAaAaAa".into(),
            "BB".into(),
            keys[0].clone(),
            keys[20].clone(),
        ]);
        let groups = group_indices(keys.iter().map(String::as_str));
        let expected = include_str!("../tests/oracle/hashmap_order.txt")
            .lines()
            .map(|s| s.parse::<usize>().unwrap())
            .collect::<Vec<_>>();
        assert_eq!(groups.iter().map(|g| g[0]).collect::<Vec<_>>(), expected);
        assert_eq!(groups.iter().map(Vec::len).sum::<usize>(), keys.len());
    }
}
