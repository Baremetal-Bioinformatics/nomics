//! Release CLI gates. NOMICS_RUNS and NOMICS_QC_BIN follow goldens.rs.
use nomics_io::bgzf::{BgzfReader, read_raw_block};
use std::fs::{self, File};
use std::io::{self, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Instant;

fn decoded(path: &Path) -> Vec<u8> {
    let mut b = Vec::new();
    BgzfReader::new(File::open(path).unwrap())
        .read_to_end(&mut b)
        .unwrap();
    b
}
fn bed_gate(ours: &Path, golden: &Path) -> Result<usize, String> {
    let a_sizes = sizes(ours);
    let b_sizes = sizes(golden);
    if a_sizes != b_sizes {
        return Err(format!(
            "BGZF block-size sequence differs: ours {} blocks, golden {}",
            a_sizes.len(),
            b_sizes.len()
        ));
    }
    // Independent decoder for the large BED streams. Production code does
    // not invoke gzip. Drain both children and check both exit statuses.
    let mut a = Command::new("gzip")
        .arg("-dc")
        .arg(ours)
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut b = Command::new("gzip")
        .arg("-dc")
        .arg(golden)
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut ar = BufReader::new(a.stdout.take().unwrap());
    let mut br = BufReader::new(b.stdout.take().unwrap());
    let mut av = [0u8; 65280];
    let mut bv = [0u8; 65280];
    let mut offset = 0;
    let mut difference = None;
    for &n in &a_sizes {
        let n = n as usize;
        if let Err(e) = ar
            .read_exact(&mut av[..n])
            .and_then(|()| br.read_exact(&mut bv[..n]))
        {
            let _ = a.kill();
            let _ = b.kill();
            let _ = a.wait();
            let _ = b.wait();
            return Err(e.to_string());
        }
        if difference.is_none()
            && let Some(i) = av[..n].iter().zip(&bv[..n]).position(|(a, b)| a != b)
        {
            difference = Some(format!(
                "decompressed byte {}: ours {}, golden {}",
                offset + i,
                av[i],
                bv[i]
            ));
        }
        offset += n;
    }
    let mut extra = [0];
    let ae = ar.read(&mut extra).unwrap();
    let be = br.read(&mut extra).unwrap();
    let aok = a.wait().unwrap().success();
    let bok = b.wait().unwrap().success();
    if !aok || !bok || ae != 0 || be != 0 {
        return Err("gzip failed or decoded size differs from BGZF ISIZE".into());
    }
    match difference {
        Some(e) => Err(e),
        None => Ok(a_sizes.len()),
    }
}
fn sizes(p: &Path) -> Vec<u32> {
    let mut r = BufReader::new(File::open(p).unwrap());
    let mut raw = Vec::new();
    let mut v = Vec::new();
    let mut at = 0;
    while read_raw_block(&mut r, &mut raw, at).unwrap() {
        v.push(u32::from_le_bytes(raw[raw.len() - 4..].try_into().unwrap()));
        at += raw.len() as u64;
    }
    v
}
fn take32(p: &[u8], at: &mut usize) -> u32 {
    let v = u32::from_le_bytes(p[*at..*at + 4].try_into().unwrap());
    *at += 4;
    v
}
fn take64(p: &[u8], at: &mut usize) -> u64 {
    let v = u64::from_le_bytes(p[*at..*at + 8].try_into().unwrap());
    *at += 8;
    v
}
fn offsets(p: &[u8]) -> Vec<(u32, u16)> {
    let mut at = 4;
    let _shift = take32(p, &mut at);
    let levels = take32(p, &mut at);
    let meta = take32(p, &mut at);
    at += meta as usize;
    let n = take32(p, &mut at);
    let mut out = Vec::new();
    let maxbin = ((1u32 << (3 * (levels + 1))) - 1) / 7;
    for _ in 0..n {
        let bins = take32(p, &mut at);
        for _ in 0..bins {
            let b = take32(p, &mut at);
            let lo = take64(p, &mut at);
            out.push((b, lo as u16));
            let chunks = take32(p, &mut at);
            for j in 0..chunks {
                let s = take64(p, &mut at);
                let e = take64(p, &mut at);
                if b < maxbin || j == 0 {
                    out.push((b, s as u16));
                    out.push((b, e as u16));
                }
            }
        }
    }
    out
}
fn same(a: &[u8], b: &[u8]) -> Result<(), String> {
    if a == b {
        Ok(())
    } else {
        let at = a
            .iter()
            .zip(b)
            .position(|(x, y)| x != y)
            .unwrap_or(a.len().min(b.len()));
        Err(format!(
            "first difference byte {at}: ours {:?}, golden {:?}; lengths {} / {}",
            a.get(at),
            b.get(at),
            a.len(),
            b.len()
        ))
    }
}
fn cli(binary: &Path, args: &[&std::ffi::OsStr]) {
    assert!(
        Command::new(binary).args(args).status().unwrap().success(),
        "CLI failed"
    );
}
#[test]
fn six_bam_release_gates() -> io::Result<()> {
    let runs = std::env::var_os("NOMICS_RUNS")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(std::env::var_os("HOME").unwrap()).join("omics/runs"));
    let binary = std::env::var_os("NOMICS_QC_BIN")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/release/nomics")
        });
    let dir = std::env::temp_dir().join(format!("nomics-mosdepth-gates-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    let bed = std::env::var_os("NOMICS_MOSDEPTH_BED")
        .map(PathBuf::from)
        .unwrap_or_else(|| runs.parent().unwrap().join("data/Agilent_v7.bed"));
    let mut failures = Vec::new();
    let mut texts = 0;
    let mut beds = 0;
    let mut csis = 0;
    for (run, sample, genome, expected) in [
        (
            "tier0-wes-20261008-r1",
            "HG002_WES",
            false,
            [50_729_515, 80_694_756],
        ),
        (
            "wgs20-20261010-r1",
            "HG002_WGS20",
            true,
            [1_091_442_954, 1_477_781_384],
        ),
        (
            "tier1-wes-20261008-r1",
            "HG002_WES",
            false,
            [5_047_709_513, 3_516_306_871],
        ),
    ] {
        if std::env::var("NOMICS_MOSDEPTH_ONLY").is_ok_and(|v| v != run) {
            continue;
        }
        let results = runs.join(run).join("results");
        let golden = results.join("reports/mosdepth").join(sample);
        for (i, (stage, kind)) in [("markduplicates", "md"), ("recalibrated", "recal")]
            .into_iter()
            .enumerate()
        {
            let name = format!("{sample}.{kind}");
            let input = results.join(format!("preprocessing/{stage}/{sample}/{name}.bam"));
            let suffixes = if genome {
                vec!["regions.bed.gz"]
            } else {
                vec!["per-base.bed.gz", "regions.bed.gz"]
            };
            let mut leaves = vec![
                "mosdepth.global.dist.txt".to_string(),
                "mosdepth.region.dist.txt".to_string(),
                "mosdepth.summary.txt".to_string(),
            ];
            for suffix in &suffixes {
                leaves.push((*suffix).to_string());
                leaves.push(format!("{suffix}.csi"));
            }
            if !input.is_file()
                || leaves
                    .iter()
                    .any(|suffix| !golden.join(format!("{name}.{suffix}")).is_file())
                || (!genome && !bed.is_file())
            {
                eprintln!("{run} {kind}: goldens absent, skipped");
                continue;
            }
            assert_eq!(
                fs::metadata(&input)?.len(),
                expected[i],
                "BAM copy incomplete: {}",
                input.display()
            );
            assert!(binary.is_file(), "build release binary first");
            let prefix = dir.join(format!("{run}-{name}"));
            let by = if genome {
                PathBuf::from("500")
            } else {
                bed.clone()
            };
            let mut cmd = Command::new(&binary);
            cmd.args(["qc", "mosdepth", "--bam"])
                .arg(&input)
                .arg("--prefix")
                .arg(&prefix)
                .arg("--by")
                .arg(by);
            if genome {
                cmd.args(["--no-per-base", "--fast-mode"]);
            }
            let start = Instant::now();
            assert!(cmd.status()?.success());
            let seconds = start.elapsed().as_secs_f64();
            eprintln!("{run} {kind}: {seconds:.3} s");
            for suffix in [
                "mosdepth.global.dist.txt",
                "mosdepth.region.dist.txt",
                "mosdepth.summary.txt",
            ] {
                texts += 1;
                let a = fs::read(prefix.with_file_name(format!(
                    "{}.{suffix}",
                    prefix.file_name().unwrap().to_string_lossy()
                )))?;
                let b = fs::read(golden.join(format!("{name}.{suffix}")))?;
                match same(&a, &b) {
                    Ok(()) => eprintln!("{run} {kind} {suffix}: byte-identical"),
                    Err(e) => failures.push(format!("{run} {kind} {suffix}: {e}")),
                }
            }
            for suffix in suffixes {
                let ours = PathBuf::from(format!("{}.{suffix}", prefix.display()));
                let g = golden.join(format!("{name}.{suffix}"));
                beds += 1;
                match bed_gate(&ours, &g) {
                    Ok(n) => eprintln!(
                        "{run} {kind} {suffix}: decompressed bytes and {n} block sizes identical"
                    ),
                    Err(e) => failures.push(format!("{run} {kind} {suffix}: {e}")),
                }
                let oc = PathBuf::from(format!("{}.csi", ours.display()));
                let gc = PathBuf::from(format!("{}.csi", g.display()));
                let rebuilt_g = dir.join("golden.csi");
                let rebuilt_o = dir.join("ours.csi");
                cli(
                    &binary,
                    &[
                        "qc".as_ref(),
                        "csi".as_ref(),
                        "--bed-gz".as_ref(),
                        g.as_os_str(),
                        "--out".as_ref(),
                        rebuilt_g.as_os_str(),
                    ],
                );
                cli(
                    &binary,
                    &[
                        "qc".as_ref(),
                        "csi".as_ref(),
                        "--bed-gz".as_ref(),
                        ours.as_os_str(),
                        "--out".as_ref(),
                        rebuilt_o.as_os_str(),
                    ],
                );
                let gb = decoded(&gc);
                let ob = decoded(&oc);
                csis += 1;
                for (part, a, b) in [
                    ("golden-side", decoded(&rebuilt_g), gb.clone()),
                    ("own-side", decoded(&rebuilt_o), ob.clone()),
                ] {
                    match same(&a, &b) {
                        Ok(()) => eprintln!(
                            "{run} {kind} {suffix}.csi {part}: decompressed payload identical"
                        ),
                        Err(e) => failures.push(format!("{run} {kind} {suffix}.csi {part}: {e}")),
                    }
                }
                match same(&fs::read(&oc)?, &fs::read(&rebuilt_o)?) {
                    Ok(()) => eprintln!(
                        "{run} {kind} {suffix}.csi own-side: complete compressed bytes identical"
                    ),
                    Err(e) => failures.push(format!(
                        "{run} {kind} {suffix}.csi own-side compressed: {e}"
                    )),
                }
                let golden_sizes = sizes(&gc) == sizes(&rebuilt_g);
                let own_sizes = sizes(&oc) == sizes(&rebuilt_o);
                eprintln!(
                    "{run} {kind} {suffix}.csi: rebuilt golden block sizes equal={golden_sizes}; rebuilt own block sizes equal={own_sizes}"
                );
                let block_sizes = sizes(&gc) == sizes(&oc);
                let within = offsets(&gb) == offsets(&ob);
                eprintln!(
                    "{run} {kind} {suffix}.csi: block-size sequences equal={block_sizes}; within-block offsets equal={within}"
                );
                // Cross-file bin compaction depends on compressed distances,
                // so these are diagnostics, not a third index normalization.
                if gb != ob {
                    eprintln!(
                        "{run} {kind} {suffix}.csi cross-file: {}",
                        same(&ob, &gb).unwrap_err()
                    );
                }
            }
        }
    }
    eprintln!(
        "mosdepth gates: {texts} text, {beds} BED, {csis} two-sided CSI; {} failures",
        failures.len()
    );
    assert!(
        failures.is_empty(),
        "{}\noutputs retained at {}",
        failures.join("\n"),
        dir.display()
    );
    fs::remove_dir_all(dir)?;
    Ok(())
}
