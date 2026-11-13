use nomics_align::ksw::{KSW_XBYTE, KSW_XSTART, KSW_XSUBO, Kswr, align2};

const MAT: [i8; 25] = [
    1, -4, -4, -4, -1, -4, 1, -4, -4, -1, -4, -4, 1, -4, -1, -4, -4, -4, 1, -1, -1, -1, -1, -1, -1,
];

fn duplicated_hit(byte_scores: bool) -> Kswr {
    // 17 bases occupy two striped vectors in either score width, leaving
    // padding lanes that must not win the qe tie or disturb score2/te2.
    let query = vec![0; 17];
    let mut target = query.clone();
    target.extend_from_slice(&[1; 40]);
    target.extend_from_slice(&query);
    align2(
        query.len() as i32,
        &query,
        target.len() as i32,
        &target,
        5,
        &MAT,
        6,
        1,
        6,
        1,
        KSW_XSUBO | KSW_XSTART | 16 | if byte_scores { KSW_XBYTE } else { 0 },
        None,
    )
}

#[test]
fn padding_lanes_preserve_score2_endpoint_and_start() {
    let expected = Kswr {
        score: 17,
        te: 16,
        qe: 16,
        score2: 17,
        te2: 73,
        tb: 0,
        qb: 0,
    };
    assert_eq!(duplicated_hit(true), expected);
    assert_eq!(duplicated_hit(false), expected);
}

#[test]
fn thread_scratch_handles_alternating_widths_and_lengths() {
    for length in [15usize, 16, 17, 31, 32, 33, 151, 17] {
        let query: Vec<u8> = (0..length).map(|i| (i % 4) as u8).collect();
        let mut target = vec![4; 11];
        target.extend_from_slice(&query);
        for byte_scores in [true, false] {
            let result = align2(
                length as i32,
                &query,
                target.len() as i32,
                &target,
                5,
                &MAT,
                6,
                1,
                6,
                1,
                KSW_XSTART | if byte_scores { KSW_XBYTE } else { 0 },
                None,
            );
            assert_eq!(result.score, length as i32);
            assert_eq!((result.qb, result.qe), (0, length as i32 - 1));
            assert_eq!((result.tb, result.te), (11, length as i32 + 10));
        }
    }
}
