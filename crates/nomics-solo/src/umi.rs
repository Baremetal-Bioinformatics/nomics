//! 1MM_All is the number of connected components of the undirected Hamming-1
//! graph, not a greedy abundance-based collapse. Read multiplicities do not
//! affect this count (SoloFeature_collapseUMI_Graph.cpp).
pub(crate) fn components(umis: &[u32]) -> u32 {
    let mut parent: Vec<usize> = (0..umis.len()).collect();
    let mut count = umis.len() as u32;
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    for (i, &u) in umis.iter().enumerate() {
        for shift in (0..24).step_by(2) {
            for delta in 1..4 {
                let other = u ^ (delta << shift);
                if other <= u {
                    continue;
                }
                if let Ok(j) = umis.binary_search(&other) {
                    let a = root(&mut parent, i);
                    let b = root(&mut parent, j);
                    if a != b {
                        parent[b] = a;
                        count -= 1;
                    }
                }
            }
        }
    }
    count
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn transitive_not_greedy() {
        // AAA -> AAC -> ACC; endpoints differ by two bases, one component.
        assert_eq!(components(&[0, 1, 5]), 1);
        assert_eq!(components(&[0, 5]), 2);
        assert_eq!(components(&[0, 1, 5, 0xffffff]), 2);
        assert_eq!(components(&[]), 0);
    }
    #[test]
    fn includes_edges_in_both_halves() {
        assert_eq!(components(&[0, 1, 1 << 22]), 1);
    }
}
