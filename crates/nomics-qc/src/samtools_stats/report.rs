use super::stats::{Stats, percentile};
use std::{fmt::Write, io};

/// C `%e` uses a signed, at-least-two-digit exponent; Rust does not.
fn scientific(x: f64) -> String {
    let s = format!("{x:.6e}");
    let Some((a, b)) = s.split_once('e') else {
        return s;
    };
    let e: i32 = b.parse().unwrap();
    format!("{a}e{e:+03}")
}
fn div32(a: u64, b: u64) -> f32 {
    if b == 0 { 0. } else { a as f32 / b as f32 }
}

pub(super) fn insert_summary(bins: &mut [[u64; 3]]) -> (usize, [u64; 3], f64, f64) {
    let mut totals = [0; 3];
    for b in bins.iter_mut() {
        for j in 0..3 {
            b[j] = (b[j] as f64 * 0.5) as u64;
            totals[j] += b[j];
        }
    }
    let mut count = totals.iter().sum::<u64>();
    let mut bulk = 0.;
    let mut avg = 0.;
    let mut end = 0;
    for (i, b) in bins.iter().enumerate() {
        let n = b.iter().sum::<u64>();
        if n > 0 {
            end = i + 1;
        }
        bulk += n as f64;
        avg += (i as u64 * n) as f64;
        if bulk / count as f64 > 0.99f32 as f64 {
            end = i + 1;
            count = bulk as u64;
            break;
        }
    }
    avg /= count.max(1) as f64;
    let mut sd = 0.;
    for (i, b) in bins.iter().enumerate().take(end).skip(1) {
        sd += b.iter().sum::<u64>() as f64 * (i as f64 - avg) * (i as f64 - avg)
            / count.max(1) as f64;
    }
    (end, totals, avg, sd.sqrt())
}
impl Stats {
    /// Finalize and format a report. Consumes the accumulator because upstream
    /// halves insert bins and mutates GC-depth bins during output.
    pub fn render(mut self, input_name: &str, reference_name: &str) -> io::Result<String> {
        self.ring.flush(-1)?;
        let (end, pairs, avg, sd) = insert_summary(&mut self.isize);
        let total = self.reads.iter().sum::<u64>();
        let avglen = div32(self.length[0], total);
        let mut s = String::new();
        writeln!(s,"# This file was produced by samtools stats (1.24+htslib-1.24) and can be plotted using plot-bamstats").unwrap();
        writeln!(s, "# This file contains statistics for all reads.").unwrap();
        writeln!(
            s,
            "# The command line was:  stats --threads 1 --reference {reference_name} {input_name}"
        )
        .unwrap();
        s.push_str("# CHK, Checksum\t[2]Read Names\t[3]Sequences\t[4]Qualities\n# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)\n");
        writeln!(
            s,
            "CHK\t{:08x}\t{:08x}\t{:08x}",
            self.chk[0], self.chk[1], self.chk[2]
        )
        .unwrap();
        s.push_str("# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.\n");
        macro_rules! sn {
            ($name:expr,$v:expr) => {
                writeln!(s, "SN\t{}:\t{}", $name, $v).unwrap()
            };
            ($name:expr,$v:expr,$comment:expr) => {
                writeln!(s, "SN\t{}:\t{}\t# {}", $name, $v, $comment).unwrap()
            };
        }
        sn!(
            "raw total sequences",
            total,
            "excluding supplementary and secondary reads"
        );
        sn!("filtered sequences", 0);
        sn!("sequences", total);
        sn!(
            "is sorted",
            u8::from(self.sorted),
            if self.sorted {
                "sorted by coordinate"
            } else {
                "not sorted by coordinate"
            }
        );
        sn!("1st fragments", self.reads[0]);
        sn!("last fragments", self.reads[1]);
        sn!("reads mapped", self.mapped);
        sn!(
            "reads mapped and paired",
            self.paired_mapped,
            "paired-end technology bit set + both mates mapped"
        );
        sn!("reads unmapped", self.unmapped);
        sn!("reads properly paired", self.proper, "proper-pair bit set");
        sn!("reads paired", self.paired, "paired-end technology bit set");
        sn!(
            "reads duplicated",
            self.dup,
            "PCR or optical duplicate bit set"
        );
        sn!("reads MQ0", self.mq0, "mapped and MQ=0");
        sn!("reads QC failed", self.qcfail);
        sn!("non-primary alignments", self.secondary);
        sn!("supplementary alignments", self.supplementary);
        sn!("total length", self.length[0], "ignores clipping");
        sn!(
            "total first fragment length",
            self.length[1],
            "ignores clipping"
        );
        sn!(
            "total last fragment length",
            self.length[2],
            "ignores clipping"
        );
        sn!("bases mapped", self.mapped_bases, "ignores clipping");
        sn!("bases mapped (cigar)", self.cigar_bases, "more accurate");
        sn!("bases trimmed", 0);
        sn!("bases duplicated", self.duplen);
        sn!("mismatches", self.nm, "from NM fields");
        sn!(
            "error rate",
            scientific(div32(self.nm, self.cigar_bases) as f64),
            "mismatches / bases mapped (cigar)"
        );
        sn!("average length", format!("{avglen:.0}"));
        sn!(
            "average first fragment length",
            format!("{:.0}", div32(self.length[1], self.reads[0]))
        );
        sn!(
            "average last fragment length",
            format!("{:.0}", div32(self.length[2], self.reads[1]))
        );
        sn!("maximum length", self.maxlen[0]);
        sn!("maximum first fragment length", self.maxlen[1]);
        sn!("maximum last fragment length", self.maxlen[2]);
        sn!(
            "average quality",
            format!(
                "{:.1}",
                if self.length[0] > 0 {
                    self.sumq / self.length[0] as f64
                } else {
                    0.
                }
            )
        );
        sn!("insert size average", format!("{avg:.1}"));
        sn!("insert size standard deviation", format!("{sd:.1}"));
        sn!("inward oriented pairs", pairs[0]);
        sn!("outward oriented pairs", pairs[1]);
        sn!("pairs with other orientation", pairs[2]);
        sn!("pairs on different chromosomes", self.anomalous / 2);
        sn!(
            "percentage of properly paired reads (%)",
            format!("{:.1}", div32(100 * self.proper, total))
        );
        let maxq = (self.maxq + usize::from(self.maxq < 255)).min(255);
        let maxlen = self.maxlen[0] + 1;
        for (g, (tag, name)) in [("FFQ", "First"), ("LFQ", "Last")].iter().enumerate() {
            writeln!(
                s,
                "# {name} Fragment Qualities. Use `grep ^{tag} | cut -f 2-` to extract this part."
            )
            .unwrap();
            s.push_str("# Columns correspond to qualities and rows to cycles. First column is the cycle number.\n");
            for i in 0..self.maxlen[g + 1] {
                write!(s, "{tag}\t{}", i + 1).unwrap();
                for q in 0..=maxq {
                    write!(s, "\t{}", self.quals[g][i][q]).unwrap();
                }
                s.push('\n');
            }
        }
        for (g, (tag, name)) in [("GCF", "first"), ("GCL", "last")].iter().enumerate() {
            writeln!(s,"# GC Content of {name} fragments. Use `grep ^{tag} | cut -f 2-` to extract this part.").unwrap();
            let mut prev = 0;
            for i in 0..200 {
                if self.gc[g][i] == self.gc[g][prev] {
                    continue;
                }
                writeln!(
                    s,
                    "{tag}\t{:.2}\t{}",
                    (i + prev) as f64 * 0.5 * 100. / 199.,
                    self.gc[g][prev]
                )
                .unwrap();
                prev = i;
            }
        }
        s.push_str("# ACGT content per cycle. Use `grep ^GCC | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N and O counts as a percentage of all A/C/G/T bases [%]\n");
        for i in 0..maxlen {
            let b = std::array::from_fn::<_, 6, _>(|j| self.bases[0][i][j] + self.bases[1][i][j]);
            base_line(&mut s, "GCC", i + 1, &b);
        }
        s.push_str("# ACGT content per cycle, read oriented. Use `grep ^GCT | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]\n");
        for i in 0..maxlen {
            base_line(&mut s, "GCT", i + 1, &self.oriented[i]);
        }
        for (g, (tag, total_tag, name)) in [("FBC", "FTC", "first"), ("LBC", "LTC", "last")]
            .iter()
            .enumerate()
        {
            writeln!(s,"# ACGT content per cycle for {name} fragments. Use `grep ^{tag} | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N and O counts as a percentage of all A/C/G/T bases [%]").unwrap();
            let mut totals = [0u64; 5];
            for i in 0..maxlen {
                for (j, t) in totals.iter_mut().enumerate() {
                    *t += self.bases[g][i][j];
                }
                base_line(&mut s, tag, i + 1, &self.bases[g][i]);
            }
            writeln!(s,"# ACGT raw counters for {name} fragments. Use `grep ^{total_tag} | cut -f 2-` to extract this part. The columns are: A,C,G,T,N base counters").unwrap();
            write!(s, "{total_tag}").unwrap();
            for t in totals {
                write!(s, "\t{t}").unwrap();
            }
            s.push('\n');
        }
        for (b, (tag, qtag)) in
            self.barcodes
                .iter()
                .zip([("BC", "QT"), ("CR", "CY"), ("OX", "BZ"), ("RX", "QX")])
        {
            if b.bases.is_empty() {
                continue;
            }
            writeln!(s,"# ACGT content per cycle for barcodes. Use `grep ^{tag}C | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N counts as a percentage of all A/C/G/T bases [%]").unwrap();
            for (i, bases) in b.bases.iter().enumerate() {
                if b.sep == Some(i) {
                    continue;
                }
                let (part, cycle) = barcode_cycle(i, b.sep);
                base_line(&mut s, &format!("{tag}C{part}"), cycle, &bases[..5]);
            }
            writeln!(
                s,
                "# Barcode Qualities. Use `grep ^{qtag}Q | cut -f 2-` to extract this part."
            )
            .unwrap();
            s.push_str("# Columns correspond to qualities and rows to barcode cycles. First column is the cycle number.\n");
            for (i, q) in b.quals.iter().enumerate() {
                if b.sep == Some(i) {
                    continue;
                }
                let (part, cycle) = barcode_cycle(i, b.sep);
                write!(s, "{qtag}Q{part}\t{cycle}").unwrap();
                if let Some(max) = b.maxq {
                    for x in &q[..=max] {
                        write!(s, "\t{x}").unwrap();
                    }
                }
                s.push('\n');
            }
        }
        s.push_str("# Insert sizes. Use `grep ^IS | cut -f 2-` to extract this part. The columns are: insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs\n");
        for (i, b) in self.isize.iter().enumerate().take(end) {
            writeln!(
                s,
                "IS\t{i}\t{}\t{}\t{}\t{}",
                b.iter().sum::<u64>(),
                b[0],
                b[1],
                b[2]
            )
            .unwrap();
        }
        for (g, (tag, title)) in [
            ("RL", "Read lengths"),
            ("FRL", "Read lengths - first fragments"),
            ("LRL", "Read lengths - last fragments"),
        ]
        .iter()
        .enumerate()
        {
            writeln!(s,"# {title}. Use `grep ^{tag} | cut -f 2-` to extract this part. The columns are: read length, count").unwrap();
            let max = if g == 0 { maxlen } else { self.maxlen[g] };
            for i in 1..=max {
                if self.lens[i][g] > 0 {
                    writeln!(s, "{tag}\t{i}\t{}", self.lens[i][g]).unwrap();
                }
            }
        }
        s.push_str("# Mapping qualities for reads !(UNMAP|SECOND|SUPPL|QCFAIL|DUP). Use `grep ^MAPQ | cut -f 2-` to extract this part. The columns are: mapq, count\n");
        for (i, &n) in self.mapq.iter().enumerate() {
            if n > 0 {
                writeln!(s, "MAPQ\t{i}\t{n}").unwrap();
            }
        }
        s.push_str("# Indel distribution. Use `grep ^ID | cut -f 2-` to extract this part. The columns are: length, number of insertions, number of deletions\n");
        for (i, b) in self.indels.iter().enumerate() {
            if b[0] > 0 || b[1] > 0 {
                writeln!(s, "ID\t{}\t{}\t{}", i + 1, b[0], b[1]).unwrap();
            }
        }
        s.push_str("# Indels per cycle. Use `grep ^IC | cut -f 2-` to extract this part. The columns are: cycle, number of insertions (fwd), .. (rev) , number of deletions (fwd), .. (rev)\n");
        for (i, b) in self.cycles.iter().enumerate() {
            if b.iter().any(|&n| n > 0) {
                writeln!(s, "IC\t{}\t{}\t{}\t{}\t{}", i + 1, b[0], b[1], b[2], b[3]).unwrap();
            }
        }
        if self.sorted {
            s.push_str(
                "# Coverage distribution. Use `grep ^COV | cut -f 2-` to extract this part.\n",
            );
            if self.ring.cov[0] > 0 {
                writeln!(s, "COV\t[<1]\t0\t{}", self.ring.cov[0]).unwrap();
            }
            for i in 1..=1000 {
                if self.ring.cov[i] > 0 {
                    writeln!(s, "COV\t[{i}-{i}]\t{i}\t{}", self.ring.cov[i]).unwrap();
                }
            }
            if self.ring.cov[1001] > 0 {
                writeln!(s, "COV\t[1000<]\t1000\t{}", self.ring.cov[1001]).unwrap();
            }
            s.push_str("# GC-depth. Use `grep ^GCD | cut -f 2-` to extract this part. The columns are: GC%, unique sequence percentiles, 10th, 25th, 50th, 75th and 90th depth percentile\n");
            let last = self.gcd.len() - 1;
            // Upstream rounds [0,igcd), sorts [0,igcd], then outputs [0,igcd).
            for b in self.gcd.iter_mut().take(last) {
                if b.1 > 0 {
                    b.0 = (100. * b.0 as f64 / b.1 as f64).round_ties_even() as f32;
                }
            }
            self.gcd
                .sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            let mut i = 0;
            while i < last {
                let gc = self.gcd[i].0;
                let mut j = i + 1;
                while j < last && (self.gcd[j].0 - gc).abs() < 0.1 {
                    j += 1;
                }
                write!(
                    s,
                    "GCD\t{gc:.1}\t{:.3}",
                    (j + 1) as f64 * 100. / self.gcd.len() as f64
                )
                .unwrap();
                for p in [10, 25, 50, 75, 90] {
                    write!(
                        s,
                        "\t{:.3}",
                        percentile(&self.gcd[i..j], p) * avglen / 20000.
                    )
                    .unwrap();
                }
                s.push('\n');
                i = j;
            }
        }
        Ok(s)
    }
}
fn barcode_cycle(i: usize, sep: Option<usize>) -> (usize, usize) {
    if let Some(s) = sep.filter(|&s| i >= s) {
        (2, i - s)
    } else {
        (1, i + 1)
    }
}
fn base_line(s: &mut String, tag: &str, cycle: usize, b: &[u64]) {
    let sum = b[..4].iter().sum::<u64>();
    if sum == 0 {
        return;
    }
    write!(s, "{tag}\t{cycle}").unwrap();
    for &n in b {
        write!(s, "\t{:.2}", 100. * n as f64 / sum as f64).unwrap();
    }
    s.push('\n');
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn c_exponent_and_rounding() {
        assert_eq!(scientific(0.002332339f32 as f64), "2.332339e-03");
        assert_eq!(scientific(0.), "0.000000e+00");
        assert_eq!(format!("{:.0}", 2.5), "2");
    }
    #[test]
    fn bulk_is_strict_float32_threshold_and_sd_omits_zero() {
        let mut b = vec![[0; 3]; 5];
        b[0] = [2, 0, 0];
        b[2] = [196, 0, 0];
        b[4] = [2, 0, 1];
        let (end, totals, avg, sd) = insert_summary(&mut b);
        // 99/100 is below the float32 representation of .99.
        assert_eq!(end, 5);
        assert_eq!(totals, [100, 0, 0]);
        assert_eq!(avg, 2.);
        assert_eq!(sd, 0.2);
        assert_eq!(b[4], [1, 0, 0]);
    }
    #[test]
    fn percentiles_use_n_plus_one_float32_interpolation() {
        let b = [(0., 10), (0., 20), (0., 40), (0., 50)];
        assert_eq!(percentile(&b, 10), 10.);
        assert_eq!(percentile(&b, 25), 12.5);
        assert_eq!(percentile(&b, 50), 30.);
        assert_eq!(percentile(&b, 90), 50.);
    }
}
