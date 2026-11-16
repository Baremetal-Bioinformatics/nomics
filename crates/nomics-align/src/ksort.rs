//! A literal safe port of the `KSORT_INIT` introsort in bwa 0.7.19's
//! `ksort.h`.  In particular, this is deliberately not Rust's stable sort:
//! callers which compare only part of an element get klib's tie permutation.

/// Port of `__ks_insertsort_##name`.
fn insertsort<T: Copy, F: FnMut(&T, &T) -> bool>(a: &mut [T], lt: &mut F) {
    for i in 1..a.len() {
        let mut j = i;
        while j > 0 && lt(&a[j], &a[j - 1]) {
            a.swap(j, j - 1);
            j -= 1;
        }
    }
}

/// Port of `ks_combsort_##name`.
fn combsort<T: Copy, F: FnMut(&T, &T) -> bool>(a: &mut [T], lt: &mut F) {
    const SHRINK_FACTOR: f64 = 1.247_330_950_103_978_7;
    let mut gap = a.len();
    loop {
        if gap > 2 {
            gap = (gap as f64 / SHRINK_FACTOR) as usize;
            if gap == 9 || gap == 10 {
                gap = 11;
            }
        }
        let mut did_swap = false;
        for i in 0..a.len() - gap {
            if lt(&a[i + gap], &a[i]) {
                a.swap(i, i + gap);
                did_swap = true;
            }
        }
        if !did_swap && gap <= 2 {
            break;
        }
    }
    if gap != 1 {
        insertsort(a, lt);
    }
}

/// Port of `ks_introsort_##name` from bwa 0.7.19 `ksort.h`.
///
/// `lt` must have the same strict-less-than meaning as the `__sort_lt` macro
/// passed to `KSORT_INIT`. `Copy` matches BWA's POD instantiations and lets
/// the port retain C's separate `rp` pivot value exactly.
pub fn introsort<T: Copy, F: FnMut(&T, &T) -> bool>(a: &mut [T], mut lt: F) {
    let n = a.len();
    if n < 1 {
        return;
    }
    if n == 2 {
        if lt(&a[1], &a[0]) {
            a.swap(0, 1);
        }
        return;
    }
    let mut d = 2_i32;
    while (1_usize << d) < n {
        d += 1;
    }
    // C allocates sizeof(size_t) * d + 2 records. It is a LIFO work stack;
    // Vec preserves the same push/pop order without raw pointers.
    let mut stack: Vec<(usize, usize, i32)> =
        Vec::with_capacity((usize::BITS as usize * d as usize) + 2);
    let mut s = 0_usize;
    let mut t = n - 1;
    d <<= 1;
    loop {
        if s < t {
            d -= 1;
            if d == 0 {
                combsort(&mut a[s..=t], &mut lt);
                t = s;
                continue;
            }
            let mut i = s;
            let mut j = t;
            let mut k = i + ((j - i) >> 1) + 1;
            if lt(&a[k], &a[i]) {
                if lt(&a[k], &a[j]) {
                    k = j;
                }
            } else {
                k = if lt(&a[j], &a[i]) { i } else { j };
            }
            let rp = a[k];
            if k != t {
                a.swap(k, t);
            }
            loop {
                i += 1;
                while lt(&a[i], &rp) {
                    i += 1;
                }
                j -= 1;
                while i <= j && lt(&rp, &a[j]) {
                    j -= 1;
                }
                if j <= i {
                    break;
                }
                a.swap(i, j);
            }
            a.swap(i, t);
            if i - s > t - i {
                if i - s > 16 {
                    stack.push((s, i - 1, d));
                }
                s = if t - i > 16 { i + 1 } else { t };
            } else {
                if t - i > 16 {
                    stack.push((i + 1, t, d));
                }
                t = if i - s > 16 { i - 1 } else { s };
            }
        } else if let Some((left, right, depth)) = stack.pop() {
            s = left;
            t = right;
            d = depth;
        } else {
            insertsort(a, &mut lt);
            return;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::introsort;
    #[test]
    fn handles_c_special_cases() {
        let mut empty: [u32; 0] = [];
        introsort(&mut empty, |a, b| a < b);
        let mut pair = [2, 1];
        introsort(&mut pair, |a, b| a < b);
        assert_eq!(pair, [1, 2]);
    }
}
