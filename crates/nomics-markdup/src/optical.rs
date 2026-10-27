//! Picard 3.4.0 optical-duplicate detection (`OpticalDuplicateFinder`,
//! `ReadNameParser`, `AbstractMarkDuplicatesCommandLineProgram.
//! trackOpticalDuplicates`). It only feeds the metrics (optical pair count,
//! library size, the optical/non-optical set histograms): with the default
//! `TAGGING_POLICY=DontTag` duplicate flags do not change.
//!
//! Not reproduced: when Picard's read-end `SortingCollection` spills to disk,
//! its codec stores x and y as 16-bit values, so coordinates above 32767 wrap
//! (`ReadEndsForMarkDuplicatesCodec` L69-70). Whether that happens depends on
//! the JVM heap. Inputs whose coordinates stay below 32768 are unaffected.

/// `OpticalDuplicateFinder.DEFAULT_OPTICAL_DUPLICATE_DISTANCE`.
pub const PIXEL_DISTANCE: i32 = 100;
/// `DEFAULT_MAX_DUPLICATE_SET_SIZE`: larger sets are not examined.
const MAX_DUPLICATE_SET_SIZE: usize = 300_000;

/// A read end's physical location (`PhysicalLocation`); `tile == -1` means
/// none (the read name did not parse).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Location {
    pub read_group: i16,
    pub tile: i16,
    pub x: i32,
    pub y: i32,
}

impl Location {
    pub const NONE: Location = Location {
        read_group: -1,
        tile: -1,
        x: -1,
        y: -1,
    };

    pub fn has_location(&self) -> bool {
        self.tile != -1
    }
}

/// `ReadNameParser.rapidParseInt`: leading digits (optional '-'); `None`
/// when there are none (Java's `NumberFormatException`).
fn rapid_parse_int(s: &[u8]) -> Option<i32> {
    let (neg, digits) = match s.first() {
        Some(b'-') => (true, &s[1..]),
        _ => (false, s),
    };
    let mut val: i32 = 0;
    let mut any = false;
    for &c in digits {
        if c.is_ascii_digit() {
            val = val.wrapping_mul(10).wrapping_add(i32::from(c - b'0'));
            any = true;
        } else {
            break;
        }
    }
    any.then_some(if neg { val.wrapping_neg() } else { val })
}

/// Tile, x and y from a read name with Picard's default (optimized) parsing:
/// the last three ':'-separated fields, accepted only when the name has 5 or
/// 7 fields (`ReadNameParser.readLocationInformation` / `getLastThreeFields`).
pub fn parse_read_name(name: &[u8]) -> Option<(i16, i32, i32)> {
    let mut tokens = [0i32; 3];
    let mut idx: i32 = 2;
    let mut fields = 0;
    let mut end = name.len();
    let mut i = name.len() as isize - 1;
    while i >= 0 && idx >= 0 {
        let at = i as usize;
        if name[at] == b':' || at == 0 {
            fields += 1;
            let start = if at == 0 { 0 } else { at + 1 };
            tokens[idx as usize] = rapid_parse_int(&name[start..end])?;
            idx -= 1;
            end = at;
        }
        i -= 1;
    }
    while i >= 0 {
        let at = i as usize;
        if name[at] == b':' || at == 0 {
            fields += 1;
        }
        i -= 1;
    }
    if fields < 3 || !(fields == 5 || fields == 7) {
        return None;
    }
    Some((tokens[0] as i16, tokens[1], tokens[2]))
}

fn close_enough(a: &Location, b: &Location, distance: i32) -> bool {
    a.has_location()
        && b.has_location()
        && a.read_group == b.read_group
        && a.tile == b.tile
        && (a.x - b.x).abs() <= distance
        && (a.y - b.y).abs() <= distance
}

fn close_enough_short(a: &Location, b: &Location, distance: i32) -> bool {
    (a.x - b.x).abs() <= distance && (a.y - b.y).abs() <= distance
}

/// `findOpticalDuplicates(list, keeper)`: `keeper` is the keeper's position
/// in `list` when it is a member (Java object identity), `keeper_given`
/// whether the caller passed a non-null keeper at all (it picks the
/// algorithm even when the keeper is not in this list).
pub fn find_optical_duplicates(
    list: &[Location],
    keeper: Option<usize>,
    keeper_given: bool,
) -> Vec<bool> {
    let n = list.len();
    let flags = vec![false; n];
    if n < 2 || n > MAX_DUPLICATE_SET_SIZE {
        return flags;
    }
    // keeperOrNull: the keeper only counts when in the list with a location.
    let actual = keeper.filter(|&k| list[k].has_location());
    if n >= if keeper_given { 4 } else { 3 } {
        flags_with_graph(list, actual, flags)
    } else {
        flags_fast(list, actual, flags)
    }
}

fn flags_fast(list: &[Location], keeper: Option<usize>, mut flags: Vec<bool>) -> Vec<bool> {
    let n = list.len();
    if let Some(k) = keeper {
        for (i, f) in flags.iter_mut().enumerate() {
            *f = i != k && close_enough(&list[k], &list[i], PIXEL_DISTANCE);
        }
    }
    for i in 0..n {
        if Some(i) == keeper {
            continue;
        }
        for j in i + 1..n {
            if Some(j) == keeper || (flags[i] && flags[j]) {
                continue;
            }
            if close_enough(&list[i], &list[j], PIXEL_DISTANCE) {
                let index = if flags[j] { i } else { j };
                flags[index] = true;
            }
        }
    }
    flags
}

fn find_rep(cluster: &mut [usize], mut node: usize) -> usize {
    let mut rep = node;
    while rep != cluster[rep] {
        rep = cluster[rep];
    }
    while node != rep {
        let next = cluster[node];
        cluster[node] = rep;
        node = next;
    }
    rep
}

fn flags_with_graph(list: &[Location], keeper: Option<usize>, mut flags: Vec<bool>) -> Vec<bool> {
    let n = list.len();
    // Reads grouped by (read group, tile); edges within a group by distance.
    let mut groups: std::collections::HashMap<i32, Vec<usize>> = std::collections::HashMap::new();
    for (i, loc) in list.iter().enumerate() {
        if loc.has_location() {
            let key = (i32::from(loc.read_group) << 16) + i32::from(loc.tile);
            groups.entry(key).or_default().push(i);
        }
    }
    let mut cluster: Vec<usize> = (0..n).collect();
    for members in groups.values() {
        for (a, &i) in members.iter().enumerate() {
            for &j in &members[a + 1..] {
                if close_enough_short(&list[i], &list[j], PIXEL_DISTANCE) {
                    let (ri, rj) = (find_rep(&mut cluster, i), find_rep(&mut cluster, j));
                    if ri != rj {
                        cluster[ri] = rj;
                    }
                }
            }
        }
    }
    let assigned: Vec<usize> = (0..n).map(|i| find_rep(&mut cluster, i)).collect();
    // Representatives per cluster; the keeper's cluster keeps the keeper.
    let mut rep: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let keeper_cluster = keeper.map(|k| assigned[k]);
    if let (Some(k), Some(c)) = (keeper, keeper_cluster) {
        rep.insert(c, k);
    }
    // Java iterates the HashMap<Integer, Integer> of node -> cluster in
    // ascending node order (small Integer keys, capacity > size).
    for i in 0..n {
        let c = assigned[i];
        match rep.get(&c).copied() {
            Some(r) if Some(i) != keeper => {
                let (cur, best) = (&list[i], &list[r]);
                if keeper_cluster != Some(c)
                    && (cur.x < best.x || (cur.x == best.x && cur.y < best.y))
                {
                    flags[r] = true;
                    rep.insert(c, i);
                } else {
                    flags[i] = true;
                }
            }
            _ => {
                rep.insert(c, i);
            }
        }
    }
    flags
}

/// Pair orientations as Picard's `ReadEnds` constants.
pub const FR: u8 = 3;
pub const RF: u8 = 5;

/// `trackOpticalDuplicates(ends, keeper, ...)` for one duplicate-pair set:
/// the number of optical duplicates. A set mixing FR and RF (by
/// `orientationForOpticalDuplicates`) is split and each half examined with
/// the same keeper. Returns `None` where Picard throws (a mixed set holding
/// another orientation).
pub fn track_optical_duplicates(
    locations: &[Location],
    orientations: &[u8],
    keeper: usize,
) -> Option<usize> {
    let has_fr = orientations.contains(&FR);
    let has_rf = orientations.contains(&RF);
    if has_fr && has_rf {
        if orientations.iter().any(|&o| o != FR && o != RF) {
            return None;
        }
        let mut total = 0;
        for want in [FR, RF] {
            let idx: Vec<usize> = (0..locations.len())
                .filter(|&i| orientations[i] == want)
                .collect();
            let sub: Vec<Location> = idx.iter().map(|&i| locations[i]).collect();
            let k = idx.iter().position(|&i| i == keeper);
            total += find_optical_duplicates(&sub, k, true)
                .iter()
                .filter(|&&f| f)
                .count();
        }
        Some(total)
    } else {
        Some(
            find_optical_duplicates(locations, Some(keeper), true)
                .iter()
                .filter(|&&f| f)
                .count(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn loc(tile: i16, x: i32, y: i32) -> Location {
        Location {
            read_group: 0,
            tile,
            x,
            y,
        }
    }

    #[test]
    fn read_names_parse_like_picard() {
        assert_eq!(
            parse_read_name(b"A00744:46:HV3GHDSXX:1:1101:1000:1200"),
            Some((1101, 1000, 1200))
        );
        assert_eq!(parse_read_name(b"M1:7:1101:15:25"), Some((1101, 15, 25)));
        assert_eq!(parse_read_name(b"SRR14724482.33"), None);
        assert_eq!(parse_read_name(b"a:b:c:1:2:3"), None); // 6 fields
        assert_eq!(parse_read_name(b"A:1:FC:1:x:2:3"), None); // tile not numeric
    }

    #[test]
    fn small_sets_use_the_keeper() {
        // Keeper and one read 50 px away: one optical duplicate.
        let l = [loc(1101, 1000, 1000), loc(1101, 1050, 1000)];
        assert_eq!(
            find_optical_duplicates(&l, Some(0), true),
            vec![false, true]
        );
        // Different tile: none.
        let l = [loc(1101, 1000, 1000), loc(1102, 1050, 1000)];
        assert_eq!(
            find_optical_duplicates(&l, Some(0), true),
            vec![false, false]
        );
    }

    #[test]
    fn graph_clusters_transitively_and_keeps_min_xy() {
        // Chain 0-1-2 within 100 px pairwise only along the chain; keeper is 3, far away.
        let l = [
            loc(1101, 1200, 1000),
            loc(1101, 1290, 1000),
            loc(1101, 1380, 1000),
            loc(1101, 5000, 5000),
        ];
        // One cluster {0,1,2}; 0 has the smallest x so 1 and 2 are optical.
        assert_eq!(
            find_optical_duplicates(&l, Some(3), true),
            vec![false, true, true, false]
        );
    }
}
