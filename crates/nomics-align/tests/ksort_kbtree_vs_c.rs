use std::path::PathBuf;
use std::process::Command;
use std::sync::OnceLock;

use nomics_align::kbtree::KBTree;
use nomics_align::ksort::introsort;

#[derive(Clone, Copy)]
struct Item {
    key: i64,
    id: u32,
    payload: [u8; 28],
}

fn next_u64(x: &mut u64) -> u64 {
    *x ^= *x << 13;
    *x ^= *x >> 7;
    *x ^= *x << 17;
    *x
}

fn oracle() -> &'static PathBuf {
    static ORACLE: OnceLock<PathBuf> = OnceLock::new();
    ORACLE.get_or_init(|| {
        let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let source = manifest.join("tests/c/harness/ksort_kbtree_oracle.c");
        let dir = std::env::temp_dir().join(format!("nomics-klib-oracle-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let binary = dir.join("ksort_kbtree_oracle");
        let status = Command::new("cc")
            .arg("-std=c11")
            .arg("-O0")
            .arg(&source)
            .arg("-o")
            .arg(&binary)
            .status()
            .unwrap();
        assert!(status.success(), "C oracle compilation failed");
        binary
    })
}

fn c_output(args: &[String]) -> String {
    let output = Command::new(oracle()).args(args).output().unwrap();
    assert!(output.status.success());
    String::from_utf8(output.stdout).unwrap()
}

fn sort_input(seed: u64, n: usize, kind: u32) -> Vec<Item> {
    let mut seed = seed;
    (0..n)
        .map(|id| {
            let r = next_u64(&mut seed);
            let key = match kind {
                0 => (r % 7) as i64 - 3,
                1 => r as i64,
                _ => {
                    if id & 1 == 1 {
                        (n / 2 + id / 2) as i64
                    } else {
                        (id / 2) as i64
                    }
                }
            };
            Item {
                key,
                id: id as u32,
                payload: [0; 28],
            }
        })
        .collect()
}

#[test]
fn introsort_permutation_matches_klib_for_ties_and_descending() {
    let mut sizes: Vec<usize> = (0..=40).collect();
    sizes.extend([127, 1024, 10_003, 100_000]);
    for &kind in &[0, 1, 2] {
        for &descending in &[false, true] {
            for &n in &sizes {
                let seed = 0x13_57_9b_df_u64 ^ n as u64 ^ ((kind as u64) << 40);
                let mut rust = sort_input(seed, n, kind);
                introsort(&mut rust, |a, b| {
                    if descending {
                        a.key > b.key
                    } else {
                        a.key < b.key
                    }
                });
                let args = [
                    "sort".into(),
                    (descending as u8).to_string(),
                    seed.to_string(),
                    n.to_string(),
                    kind.to_string(),
                ];
                let c_ids: Vec<u32> = c_output(&args)
                    .split_whitespace()
                    .map(str::parse)
                    .collect::<Result<_, _>>()
                    .unwrap();
                assert_eq!(
                    rust.iter().map(|x| x.id).collect::<Vec<_>>(),
                    c_ids,
                    "kind={kind} desc={descending} n={n}"
                );
            }
        }
    }
}

#[test]
fn c_instrumentation_confirms_combsort_depth_fallback() {
    let out = c_output(&["fallback".into(), "4096".into()]);
    assert!(
        out.trim().parse::<u32>().unwrap() > 0,
        "median-of-three killer did not reach ks_combsort"
    );
}

fn cmp(a: &Item, b: &Item) -> i32 {
    (a.key > b.key) as i32 - (a.key < b.key) as i32
}

#[test]
fn btree_duplicate_queries_mutation_and_traversal_match_klib() {
    for &size in &[512_usize, 200] {
        // 40-byte C key gives t=5 and t=2.
        for case in 0..24_u64 {
            let seed0 = 0x9e37_79b9_7f4a_7c15_u64 ^ case;
            let ops = 150 + (case as usize % 17);
            let mut tree = KBTree::new(size, 40, cmp).unwrap();
            assert_eq!(tree.degree(), ((size - 12) / 48 + 1) >> 1);
            let mut seed = seed0;
            let mut lines = String::new();
            for id in 0..ops {
                let key = (next_u64(&mut seed) % 23) as i64 - 11;
                let query = Item {
                    key,
                    id: 0,
                    payload: [0; 28],
                };
                let (lower, upper) = tree.interval_locations(&query);
                let got = tree.getp(&query).map(|x| x.id);
                assert_eq!(
                    got,
                    lower
                        .and_then(|_| tree.intervalp(&query).0)
                        .filter(|x| x.key == key)
                        .map(|x| x.id)
                );
                if lower.is_some() {
                    assert_eq!(
                        tree.interval_lower_mut(&query).map(|x| x.id),
                        tree.intervalp(&query).0.map(|x| x.id)
                    );
                }
                let lower_id = lower
                    .and_then(|_| tree.intervalp(&query).0.map(|x| x.id))
                    .map_or(-1, |x| x as i32);
                let upper_id = upper
                    .and_then(|_| tree.intervalp(&query).1.map(|x| x.id))
                    .map_or(-1, |x| x as i32);
                lines.push_str(&format!("Q {lower_id} {upper_id}\n"));
                if let Some(p) = lower {
                    if (tree.get_location_mut(p).id ^ next_u64(&mut seed) as u32) & 1 != 0 {
                        tree.get_location_mut(p).payload[0] ^= 0x5a;
                    } else {
                        tree.putp(Item {
                            key,
                            id: id as u32,
                            payload: [0; 28],
                        });
                    }
                } else {
                    tree.putp(Item {
                        key,
                        id: id as u32,
                        payload: [0; 28],
                    });
                }
            }
            lines.push_str("T ");
            tree.traverse(|x| lines.push_str(&format!("{} ", x.id)));
            lines.push('\n');
            let expected = c_output(&[
                "btree".into(),
                size.to_string(),
                seed0.to_string(),
                ops.to_string(),
            ]);
            assert_eq!(lines, expected, "size={size} case={case}");
        }
    }
}
