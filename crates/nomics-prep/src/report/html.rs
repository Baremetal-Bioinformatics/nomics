//! Literal HTML/log layout independently implemented from fastp v1.1.0.
use super::{ReadStats, Report, float, list};
use std::fmt::Write as _;

fn row(s: &mut String, key: &str, value: impl std::fmt::Display) {
    writeln!(
        s,
        "<tr><td class='col1'>{key}</td><td class='col2'>{value}</td></tr>"
    )
    .unwrap();
}

fn number(n: u64) -> String {
    let mut v = n as f64;
    let units = ["", "K", "M", "G", "T", "P"];
    let mut order = 0;
    while v > 1000.0 {
        order += 1;
        v /= 1000.0;
    }
    if order == 0 {
        n.to_string()
    } else {
        format!("{v:.6} {}", units[order])
    }
}

fn fixed(v: f64) -> String {
    if v.is_finite() {
        format!("{v:.6}")
    } else {
        float(v)
    }
}

fn percent(n: u64, d: u64) -> String {
    if d == 0 {
        "0.0".into()
    } else {
        format!("{:.6}", n as f64 * 100.0 / d as f64)
    }
}

fn section(s: &mut String, id: &str, title: &str, class: &str) {
    writeln!(
        s,
        "<div class='{class}' onclick=showOrHide('{id}')>{title}</div>\n<div id='{id}'>"
    )
    .unwrap();
}

fn kmer(mut value: usize, len: usize) -> String {
    let mut bases = vec![b'A'; len];
    for b in bases.iter_mut().rev() {
        *b = b"ATCG"[value & 3];
        value >>= 2;
    }
    String::from_utf8(bases).unwrap()
}

impl ReadStats {
    fn html_kmers(&self, s: &mut String, stage: &str, mate: &str) {
        let title = format!("{stage}: {mate}: KMER counting");
        let id = title.replace([' ', ':'], "_");
        writeln!(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('{id}')>{title}</a></div>\n<div  id='{id}'>\n<div class='sub_section_tips'>Darker background means larger counts. The count will be shown on mouse over.</div>\n<table class='kmer_table' style='width:680px;'>").unwrap();
        s.push_str("<tr><td></td>");
        for j in 0..16 {
            write!(s, "<td style='color:#333333'>{}</td>", kmer(j, 2)).unwrap();
        }
        s.push_str("</tr>\n");
        for i in 0..64 {
            write!(s, "<tr><td style='color:#333333'>{}</td>", kmer(i, 3)).unwrap();
            for j in 0..16 {
                let idx = i * 16 + j;
                let count = self.kmers[idx];
                // fastp allocates 2048 slots but prints only the first 1024.
                let prop = count as f64 / ((self.bases + 1) as f64 / 2048.0);
                let frac = if prop > 2.0 {
                    (prop - 2.0) / 20.0 + 0.5
                } else if prop < 0.5 {
                    prop
                } else {
                    0.5
                };
                let color = ((1.0 - frac.clamp(0.01, 1.0)) * 255.0) as u32;
                let word = kmer(idx, 5);
                write!(s, "<td style='background:#{color:02x}{color:02x}{color:02x}' title='{word}: {count}\n{} times as mean value'>{word}</td>", float(prop)).unwrap();
            }
            s.push_str("</tr>\n");
        }
        s.push_str("</table>\n</div>\n");
    }

    fn html_curve(&self, s: &mut String, stage: &str, mate: &str, quality: bool) {
        let kind = if quality { "quality" } else { "base contents" };
        let title = format!("{stage}: {mate}: {kind}");
        let id = title.replace([' ', ':'], "_");
        writeln!(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('{id}')>{title}</a></div>\n<div id='{id}'>\n<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n<div class='figure' id='plot_{id}'></div>\n</div>\n\n<script type=\"text/javascript\">").unwrap();
        let cycles = self.cycles.len();
        let mut x: Vec<usize> = (1..=cycles.min(if cycles > 300 { 40 } else { cycles })).collect();
        if cycles > 300 {
            let mut pos = 40.0;
            loop {
                pos *= 1.05;
                if pos >= cycles as f64 {
                    break;
                }
                x.push(pos as usize);
            }
            if x.last().copied() != Some(cycles) {
                x.push(cycles);
            }
        }
        let names = if quality {
            &["A", "T", "C", "G", "mean"][..]
        } else {
            &["A", "T", "C", "G", "N", "GC"][..]
        };
        let mut colors = vec![
            "128,128,0,1.0",
            "128,0,128,1.0",
            "0,255,0,1.0",
            "0,0,255,1.0",
        ];
        if !quality {
            colors.push("255, 0, 0, 1.0");
        }
        colors.push("20,20,20,1.0");
        s.push_str("var data=[");
        for (base, color) in names.iter().zip(colors) {
            s.push_str("{x:[");
            list(s, x.iter().map(usize::to_string));
            s.push_str("],y:[");
            let curve = self.curve(base, quality);
            let mut previous = 0;
            list(
                s,
                x.iter().map(|&end| {
                    let sum: f64 = curve[previous..end].iter().sum();
                    let n = end - previous;
                    previous = end;
                    if n == 0 {
                        "0.0".into()
                    } else {
                        float(sum / n as f64)
                    }
                }),
            );
            let name = if quality {
                base.to_string()
            } else {
                let count: u64 = self
                    .cycles
                    .iter()
                    .map(|c| {
                        if *base == "GC" {
                            c.count[3] + c.count[7]
                        } else {
                            c.count[(base.as_bytes()[0] & 7) as usize]
                        }
                    })
                    .sum();
                let mut p = if self.bases == 0 {
                    "-nan".into()
                } else {
                    percent(count, self.bases)
                };
                p.truncate(5);
                format!("{base}({p}%)")
            };
            write!(
                s,
                "],name: '{name}',mode:'lines',line:{{color:'rgba({color})', width:1}}\n}},"
            )
            .unwrap();
        }
        writeln!(s, "];\nvar layout={{title:'', xaxis:{{title:'position'{}}}, yaxis:{{title:'{}'}}}};\nPlotly.newPlot('plot_{id}', data, layout);\n</script>", if cycles > 300 { ",type:'log'" } else { "" }, if quality { "quality" } else { "base content ratios" }).unwrap();
    }
}

impl Report {
    /// Literal fastp HTML. Timestamps are explicit so tests can reproduce a
    /// reference's two independently sampled title/footer times.
    pub fn html(
        &self,
        stats: &crate::Stats,
        command: &str,
        title_time: &str,
        footer_time: &str,
    ) -> String {
        let mut s = format!(
            "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" /><title>fastp report at {title_time} </title>{}</head><body><div id='container'>\n<h1 style='text-align:left;'><a href='https://github.com/OpenGene/fastp' target='_blank' style='color:#663355;text-decoration:none;'>fastp report</a>\n<div class='section_div'>\n",
            include_str!("head.html")
        );
        section(
            &mut s,
            "summary",
            "<a name='summary'>Summary</a>",
            "section_title",
        );
        section(&mut s, "general", "General", "subsection_title");
        s.push_str("<table class='summary_table'>\n");
        row(
            &mut s,
            "fastp version:",
            "1.1.0 (<a href='https://github.com/OpenGene/fastp'>https://github.com/OpenGene/fastp</a>)",
        );
        row(
            &mut s,
            "sequencing:",
            format!(
                "paired end ({} cycles + {} cycles)",
                self.reads[0].cycles.len(),
                self.reads[1].cycles.len()
            ),
        );
        for (part, label) in ["before", "after"].into_iter().enumerate() {
            row(
                &mut s,
                &format!("mean length {label} filtering:"),
                format!(
                    "{}bp, {}bp",
                    self.reads[part * 2].mean_length(),
                    self.reads[part * 2 + 1].mean_length()
                ),
            );
        }
        row(
            &mut s,
            "duplication rate:",
            format!("{:.6}%", self.duplicate.rate() * 100.0),
        );
        row(&mut s, "Insert size peak:", self.peak());
        s.push_str("</table>\n</div>\n");
        for (part, label) in ["Before", "After"].into_iter().enumerate() {
            section(
                &mut s,
                &format!("{}_filtering_summary", label.to_lowercase()),
                &format!("{label} filtering"),
                "subsection_title",
            );
            s.push_str("<table class='summary_table'>\n");
            let a = &self.reads[part * 2];
            let b = &self.reads[part * 2 + 1];
            let bases = a.bases + b.bases;
            row(&mut s, "total reads:", number(a.reads + b.reads));
            row(&mut s, "total bases:", number(bases));
            for (i, q) in [20, 30, 40].into_iter().enumerate() {
                row(
                    &mut s,
                    &format!("Q{q} bases:"),
                    format!(
                        "{} ({}%)",
                        number(a.q[i] + b.q[i]),
                        percent(a.q[i] + b.q[i], bases)
                    ),
                );
            }
            row(
                &mut s,
                "GC content:",
                format!("{}%", percent(a.gc + b.gc, bases)),
            );
            s.push_str("</table>\n</div>\n");
        }
        section(
            &mut s,
            "filtering_result",
            "Filtering result",
            "subsection_title",
        );
        s.push_str("<table class='summary_table'>\n");
        for (label, count) in [
            ("reads passed filters:", stats.passed_reads),
            ("reads with low quality:", stats.low_quality_reads),
            ("reads with too many N:", stats.too_many_n_reads),
            ("reads too short:", stats.too_short_reads),
        ] {
            row(
                &mut s,
                label,
                format!(
                    "{} ({}%)",
                    number(count),
                    fixed(count as f64 * 100.0 / (stats.pairs_in * 2) as f64)
                ),
            );
        }
        s.push_str("</table>\n</div>\n</div>\n</div>\n<div class='section_div'>\n");
        section(
            &mut s,
            "insert_size",
            "<a name='summary'>Insert size estimation</a>",
            "section_title",
        );
        let limit =
            (self.reads[0].cycles.len() + self.reads[1].cycles.len()).saturating_sub(30) as i64;
        let total = limit.clamp(1, 512) as usize;
        let all: f64 =
            self.insert[..total].iter().map(|&n| n as f64).sum::<f64>() + self.insert[512] as f64;
        let unknown = fixed(self.insert[512] as f64 * 100.0 / all);
        write!(s, "<div id='insert_size_figure'>\n<div class='figure' id='plot_insert_size' style='height:400px;'></div>\n</div>\n<div class='sub_section_tips'>This estimation is based on paired-end overlap analysis, and there are {unknown}% reads found not overlapped. <br /> The nonoverlapped read pairs may have insert size &lt;30 or &gt;{}, or contain too much sequencing errors to be detected as overlapped.</div>\n\n<script type=\"text/javascript\">\nvar data=[{{x:[", limit.max(1)).unwrap();
        list(&mut s, (0..total).map(|n| n.to_string()));
        s.push_str("],y:[");
        list(
            &mut s,
            self.insert[..total].iter().map(|&n| {
                float(if all > 0.0 {
                    n as f64 * 100.0 / all
                } else {
                    0.0
                })
            }),
        );
        writeln!(s, "],name: 'Percent (%)  ',type:'bar',line:{{color:'rgba(128,0,128,1.0)', width:1}}\n}}];\nvar layout={{title:'Insert size distribution ({unknown}% reads are with unknown length)', xaxis:{{title:'Insert size'}}, yaxis:{{title:'Read percent (%)'}}}};\nPlotly.newPlot('plot_insert_size', data, layout);\n</script>").unwrap();
        s.push_str("</div>\n</div>\n<div class='section_div'>\n<div class='section_title' onclick=showOrHide('quality_stat')><a name='summary'>Filtering statistics</a></div>\n<table id='quality_stat' class='section_table'>\n");
        for quality in [true, false] {
            if !quality {
                s.push_str("<tr style='height:20px;background:#999999;'></tr>\n");
            }
            for mate in 0..2 {
                s.push_str("<tr><td>\n");
                self.reads[mate].html_curve(
                    &mut s,
                    "Before filtering",
                    &format!("read{}", mate + 1),
                    quality,
                );
                s.push_str("</td><td>\n");
                self.reads[mate + 2].html_curve(
                    &mut s,
                    "After filtering",
                    &format!("read{}", mate + 1),
                    quality,
                );
                s.push_str("</td></tr>\n");
            }
            if quality {
                s.push_str("<tr style='height:20px;background:#999999;'></tr>\n<tr><td>\n");
                self.html_hist(&mut s, 0);
                s.push_str("</td><td>\n");
                self.html_hist(&mut s, 1);
                s.push_str("</td></tr>\n");
            }
        }
        s.push_str("<tr style='height:20px;background:#999999;'></tr>\n");
        for mate in 0..2 {
            s.push_str("<tr><td>\n");
            self.reads[mate].html_kmers(&mut s, "Before filtering", &format!("read{}", mate + 1));
            s.push_str("</td><td>\n");
            self.reads[mate + 2].html_kmers(
                &mut s,
                "After filtering",
                &format!("read{}", mate + 1),
            );
            s.push_str("</td></tr>\n");
        }
        write!(s, "</table>\n</div>\n\n</div>\n<div id='footer'> <p>{command}</p>fastp 1.1.0, at {footer_time} </div></body></html>").unwrap();
        s
    }

    fn html_hist(&self, s: &mut String, part: usize) {
        let title = format!(
            "{} filtering: quality score histogram",
            if part == 0 { "Before" } else { "After" }
        );
        let id = title.replace([' ', ':'], "_");
        writeln!(s, "<div class='figure' style='height:400' id='plot_{id}'></div>\n<script language='javascript'>").unwrap();
        for mate in 0..2 {
            write!(s, " var hist{} = {{\n x: [", mate + 1).unwrap();
            let hist = &self.reads[part * 2 + mate].qual_hist;
            list(
                s,
                (33..128)
                    .filter(|&q| hist[q] > 0)
                    .map(|q| (q - 33).to_string()),
            );
            s.push_str("],\n y: [");
            list(s, hist[33..].iter().filter(|&&n| n > 0).map(u64::to_string));
            writeln!(s, "],\nname: 'read{}',\ntype: 'bar',\n}};", mate + 1).unwrap();
        }
        writeln!(s, "var data = [hist1, hist2];\nvar layout = {{barmode: 'group', title: '{title}',xaxis:{{title:'base quality score'}}, yaxis:{{title:'base count'}}}}; \nPlotly.newPlot('plot_{id}', data, layout); \n</script>").unwrap();
    }

    /// fastp's stderr summary, with report-path provenance supplied by caller.
    pub fn log(
        &self,
        stats: &crate::Stats,
        command: &str,
        json_path: &str,
        html_path: &str,
        seconds: u64,
    ) -> String {
        let mut s = String::new();
        for (i, label) in [
            "Read1 before filtering",
            "Read2 before filtering",
            "Read1 after filtering",
            "Read2 after filtering",
        ]
        .into_iter()
        .enumerate()
        {
            let r = &self.reads[i];
            writeln!(
                s,
                "{label}:\ntotal reads: {}\ntotal bases: {}",
                r.reads, r.bases
            )
            .unwrap();
            for (q, n) in [20, 30, 40].into_iter().zip(r.q) {
                writeln!(
                    s,
                    "Q{q} bases: {n}({}%)",
                    float(n as f64 * 100.0 / r.bases as f64)
                )
                .unwrap();
            }
            s.push('\n');
        }
        write!(s, "Filtering result:\nreads passed filter: {}\nreads failed due to low quality: {}\nreads failed due to too many N: {}\nreads failed due to too short: {}\n\nDuplication rate: {}%\n\nInsert size peak (evaluated by paired-end reads): {}\n\nJSON report: {json_path}\nHTML report: {html_path}\n\n{command}\nfastp v1.1.0, time used: {seconds} seconds\n",stats.passed_reads, stats.low_quality_reads, stats.too_many_n_reads, stats.too_short_reads, float(self.duplicate.rate()*100.0), self.peak()).unwrap();
        s
    }
}

/// UTC wall-clock timestamp in fastp's literal spelling. fastp uses its process
/// local timezone; only the two timestamp values are proposed gate exemptions.
pub fn timestamp_now() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    timestamp(secs)
}

fn timestamp(secs: u64) -> String {
    let mut days = secs / 86400;
    let mut year = 1970u64;
    let leap = |y: u64| y.is_multiple_of(4) && (!y.is_multiple_of(100) || y.is_multiple_of(400));
    loop {
        let n = if leap(year) { 366 } else { 365 };
        if days < n {
            break;
        }
        days -= n;
        year += 1;
    }
    let months = [
        31,
        if leap(year) { 29 } else { 28 },
        31,
        30,
        31,
        30,
        31,
        31,
        30,
        31,
        30,
        31,
    ];
    let mut month = 1;
    for n in months {
        if days < n {
            break;
        }
        days -= n;
        month += 1;
    }
    format!(
        "{year}-{month:02}-{:02}      {:02}:{:02}:{:02}",
        days + 1,
        secs / 3600 % 24,
        secs / 60 % 60,
        secs % 60
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn utc_dates() {
        assert_eq!(timestamp(0), "1970-01-01      00:00:00");
        assert_eq!(timestamp(951868799), "2000-02-29      23:59:59");
        assert_eq!(timestamp(4107542400), "2100-03-01      00:00:00");
    }
}
