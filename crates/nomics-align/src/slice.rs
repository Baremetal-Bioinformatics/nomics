//! Bounded-memory paired FASTQ selection using the exact BWA-MEM batch path.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Cursor, Read as IoRead, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::mpsc::sync_channel;
use std::time::Instant;

use nomics_io::bam::RecordView;
use nomics_io::bgzf::{BgzfWriter, Profile};
use nomics_io::deflate::HtsjdkDeflater;

use crate::index::Index;
use crate::mem::{self, GpuSeeds};
use crate::sam::{BamRecords, Read};

pub struct Args {
    pub index: PathBuf,
    pub fastq1: PathBuf,
    pub fastq2: PathBuf,
    pub out1: PathBuf,
    pub out2: PathBuf,
    pub contig: String,
    pub dense_sa_cache: PathBuf,
    pub dense_sa_interval: u64,
    pub threads: usize,
}

struct Input {
    reader: BufReader<Box<dyn IoRead + Send>>,
    child: Option<Child>,
}

impl Input {
    fn open(path: &Path) -> io::Result<Self> {
        if path.extension().is_some_and(|ext| ext == "gz") {
            // GzReader currently reads the whole compressed file. `gzip -dc`
            // supplies an incremental stream for WGS-sized inputs instead.
            let mut child = Command::new("gzip")
                .arg("-dc")
                .arg("--")
                .arg(path)
                .stdout(Stdio::piped())
                .spawn()?;
            let stdout = child.stdout.take().unwrap();
            Ok(Self {
                reader: BufReader::with_capacity(1 << 20, Box::new(stdout)),
                child: Some(child),
            })
        } else {
            Ok(Self {
                reader: BufReader::with_capacity(1 << 20, Box::new(File::open(path)?)),
                child: None,
            })
        }
    }

    fn finish(mut self) -> io::Result<()> {
        if let Some(mut child) = self.child.take() {
            let status = child.wait()?;
            if !status.success() {
                return Err(io::Error::other(format!("gzip exited with {status}")));
            }
        }
        Ok(())
    }
}

enum Output {
    Plain(BufWriter<File>),
    Gzip(BgzfWriter<BufWriter<File>, HtsjdkDeflater>),
}

impl Output {
    fn open(path: &Path) -> io::Result<Self> {
        let file = BufWriter::with_capacity(1 << 20, File::create(path)?);
        if path.extension().is_some_and(|ext| ext == "gz") {
            Ok(Self::Gzip(BgzfWriter::new(
                file,
                HtsjdkDeflater::new(),
                Profile::HTSJDK,
            )))
        } else {
            Ok(Self::Plain(file))
        }
    }

    fn write_all(&mut self, data: &[u8]) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.write_all(data),
            Self::Gzip(w) => w.write_all(data),
        }
    }

    fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(mut w) => w.flush(),
            Self::Gzip(w) => w.finish()?.flush(),
        }
    }
}

struct Batch {
    first_pair: u64,
    reads: Vec<Read>,
    raw1: Vec<u8>,
    raw2: Vec<u8>,
    ends1: Vec<usize>,
    ends2: Vec<usize>,
    seeds: Option<GpuSeeds>,
}

fn raw_read(reader: &mut impl BufRead, raw: &mut Vec<u8>) -> io::Result<Option<Read>> {
    let start = raw.len();
    for line in 0..4 {
        if reader.read_until(b'\n', raw)? == 0 {
            raw.truncate(start);
            return if line == 0 {
                Ok(None)
            } else {
                Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "FASTQ record ends before its fourth line",
                ))
            };
        }
    }
    mem::read_fastq(&mut Cursor::new(&raw[start..]), &mut Vec::new())?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "empty FASTQ record"))
        .map(Some)
}

fn read_batch(r1: &mut impl BufRead, r2: &mut impl BufRead, first_pair: u64) -> io::Result<Batch> {
    let mut batch = Batch {
        first_pair,
        reads: Vec::new(),
        raw1: Vec::new(),
        raw2: Vec::new(),
        ends1: Vec::new(),
        ends2: Vec::new(),
        seeds: None,
    };
    let mut bases = 0i64;
    loop {
        let a = raw_read(r1, &mut batch.raw1)?;
        let b = raw_read(r2, &mut batch.raw2)?;
        match (a, b) {
            (Some(a), Some(b)) => {
                if a.name != b.name {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "FASTQ mate names differ at pair {}",
                            first_pair + batch.ends1.len() as u64
                        ),
                    ));
                }
                bases += (a.seq.len() + b.seq.len()) as i64;
                batch.reads.extend([a, b]);
                batch.ends1.push(batch.raw1.len());
                batch.ends2.push(batch.raw2.len());
                if bases >= mem::CHUNK_SIZE {
                    break;
                }
            }
            (None, None) => break,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "R1 and R2 differ in read count",
                ));
            }
        }
    }
    Ok(batch)
}

#[derive(Default)]
struct Stats {
    pairs: u64,
    primary: u64,
    other: u64,
    batches: u64,
}

/// Records are emitted in pair order. Each pair begins with the primary R1
/// record, including when it is unmapped; this also handles repeated names.
fn select(
    batch: &Batch,
    records: &mut BamRecords,
    target: i32,
) -> io::Result<(Vec<u8>, Vec<u8>, u64, u64)> {
    if let Some(error) = records.take_error() {
        return Err(error);
    }
    let pairs = batch.reads.len() / 2;
    let mut hits = vec![0u8; pairs];
    let mut pair = None;
    let mut at = 0;
    while at < records.bytes.len() {
        let len =
            u32::from_le_bytes(records.bytes[at..at + 4].try_into().map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "truncated BAM record size")
            })?) as usize;
        let end = at + 4 + len;
        let raw = records
            .bytes
            .get(at + 4..end)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "truncated BAM record"))?;
        let record = RecordView::new(raw)?;
        let flag = record.flag();
        if flag & 0x40 != 0 && flag & 0x900 == 0 {
            let next = pair.map_or(0, |p| p + 1);
            if next >= pairs || record.read_name() != batch.reads[2 * next].name {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "BAM pair order/name differs from FASTQ",
                ));
            }
            pair = Some(next);
        }
        let p = pair.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "BAM pair has no primary R1")
        })?;
        if record.read_name() != batch.reads[2 * p].name {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "BAM read name differs from FASTQ",
            ));
        }
        if record.ref_id() == target || record.mate_ref_id() == target {
            let primary = flag & 0x900 == 0 && flag & 4 == 0 && record.ref_id() == target;
            hits[p] = hits[p].max(if primary { 2 } else { 1 });
        }
        at = end;
    }
    if pair.map_or(0, |p| p + 1) != pairs {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "BAM has fewer pairs than FASTQ",
        ));
    }
    let mut out1 = Vec::new();
    let mut out2 = Vec::new();
    let (mut primary, mut other) = (0, 0);
    for (p, hit) in hits.into_iter().enumerate() {
        if hit == 0 {
            continue;
        }
        primary += u64::from(hit == 2);
        other += u64::from(hit == 1);
        let s1 = if p == 0 { 0 } else { batch.ends1[p - 1] };
        let s2 = if p == 0 { 0 } else { batch.ends2[p - 1] };
        out1.extend_from_slice(&batch.raw1[s1..batch.ends1[p]]);
        out2.extend_from_slice(&batch.raw2[s2..batch.ends2[p]]);
    }
    Ok((out1, out2, primary, other))
}

pub fn run(args: &Args) -> io::Result<()> {
    if args.threads == 0
        || args.fastq1 == args.out1
        || args.fastq1 == args.out2
        || args.fastq2 == args.out1
        || args.fastq2 == args.out2
        || args.out1 == args.out2
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "threads must be positive and input/output paths distinct",
        ));
    }
    let mut idx = Index::load(&args.index)?;
    idx.attach_dense_sa_with_interval(&args.dense_sa_cache, args.dense_sa_interval, args.threads)?;
    let target = idx
        .bns
        .anns
        .iter()
        .position(|ann| ann.name == args.contig)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("contig {} absent from index", args.contig),
            )
        })? as i32;
    let mut input1 = Input::open(&args.fastq1)?;
    let mut input2 = Input::open(&args.fastq2)?;
    let output1 = Output::open(&args.out1)?;
    let output2 = Output::open(&args.out2)?;
    let gpu = std::env::var_os("NOMICS_GPU_SEED").is_some();
    let inflight = if args.threads > 1 { 2 } else { 1 };
    let workers = (args.threads / inflight).max(1);
    let started = Instant::now();
    let idx_ref = &idx;
    let stats = std::thread::scope(|scope| -> io::Result<Stats> {
        let (batch_tx, batch_rx) = sync_channel::<Batch>(inflight);
        let producer = scope.spawn(move || -> io::Result<()> {
            let mut first_pair = 0u64;
            loop {
                let mut batch = read_batch(&mut input1.reader, &mut input2.reader, first_pair)?;
                if batch.reads.is_empty() {
                    break;
                }
                first_pair += (batch.reads.len() / 2) as u64;
                batch.seeds = gpu.then(|| mem::try_gpu_seeds(idx_ref, &batch.reads)).flatten();
                if batch_tx.send(batch).is_err() {
                    return Ok(());
                }
            }
            input1.finish()?;
            input2.finish()
        });
        let (write_tx, write_rx) = sync_channel::<(Vec<u8>, Vec<u8>)>(2);
        let writer = scope.spawn(move || -> io::Result<()> {
            let (mut out1, mut out2) = (output1, output2);
            for (a, b) in write_rx {
                out1.write_all(&a)?;
                out2.write_all(&b)?;
            }
            out1.finish()?;
            out2.finish()
        });
        let mut stats = Stats::default();
        let mut pending = VecDeque::new();
        let mut result = Ok(());
        for mut batch in &batch_rx {
            pending.push_back(scope.spawn(move || {
                let mut records = BamRecords::default();
                mem::process_batch_bam(
                    idx_ref,
                    b"",
                    &batch.reads,
                    batch.seeds.take(),
                    batch.first_pair,
                    workers,
                    &mut records,
                );
                (batch, records)
            }));
            if pending.len() == inflight {
                result = finish_batch(
                    pending.pop_front().unwrap(),
                    target,
                    &write_tx,
                    &mut stats,
                    started,
                );
                if result.is_err() {
                    break;
                }
            }
        }
        drop(batch_rx);
        if result.is_ok() {
            while let Some(handle) = pending.pop_front() {
                result = finish_batch(handle, target, &write_tx, &mut stats, started);
                if result.is_err() {
                    break;
                }
            }
        }
        drop(write_tx);
        let producer_result = producer
            .join()
            .map_err(|_| io::Error::other("slice reader panicked"))?;
        let writer_result = writer
            .join()
            .map_err(|_| io::Error::other("slice writer panicked"))?;
        result?;
        producer_result?;
        writer_result?;
        Ok(stats)
    })?;
    eprintln!(
        "slice total: processed={} kept={} primary={} secondary/supplementary/mate_only={} elapsed={:.1}s rate={:.0} pairs/s",
        stats.pairs,
        stats.primary + stats.other,
        stats.primary,
        stats.other,
        started.elapsed().as_secs_f64(),
        stats.pairs as f64 / started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn finish_batch(
    handle: std::thread::ScopedJoinHandle<'_, (Batch, BamRecords)>,
    target: i32,
    tx: &std::sync::mpsc::SyncSender<(Vec<u8>, Vec<u8>)>,
    stats: &mut Stats,
    started: Instant,
) -> io::Result<()> {
    let (batch, mut records) = handle
        .join()
        .map_err(|_| io::Error::other("slice alignment panicked"))?;
    if batch.first_pair != stats.pairs {
        return Err(io::Error::other("slice batches completed out of order"));
    }
    let (a, b, primary, other) = select(&batch, &mut records, target)?;
    tx.send((a, b))
        .map_err(|_| io::Error::other("slice writer stopped"))?;
    stats.pairs += (batch.reads.len() / 2) as u64;
    stats.primary += primary;
    stats.other += other;
    stats.batches += 1;
    if stats.batches % 10 == 0 {
        eprintln!(
            "slice: processed={} kept={} rate={:.0} pairs/s",
            stats.pairs,
            stats.primary + stats.other,
            stats.pairs as f64 / started.elapsed().as_secs_f64()
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    fn add_record(records: &mut BamRecords, name: &[u8], flag: u16, rid: i32, mate_rid: i32) {
        let mut raw = vec![0u8; 32 + name.len() + 1];
        raw[..4].copy_from_slice(&rid.to_le_bytes());
        raw[8] = (name.len() + 1) as u8;
        raw[14..16].copy_from_slice(&flag.to_le_bytes());
        raw[20..24].copy_from_slice(&mate_rid.to_le_bytes());
        raw[32..32 + name.len()].copy_from_slice(name);
        records
            .bytes
            .extend_from_slice(&(raw.len() as u32).to_le_bytes());
        records.bytes.extend_from_slice(&raw);
    }

    #[test]
    fn selects_all_record_classes_and_preserves_fastq_bytes() {
        let r1 = b"@one/1 comment a\r\nACGT\r\n+one\r\n!!!!\r\n@two/1 comment b\nTT\n+\n##\n@three/1 mate\nCC\n+\n%%\n";
        let r2 = b"@one/2 comment c\nTGCA\n+\n!!!!\n@two/2 comment d\nAA\n+\n$$\n@three/2 mate\nGG\n+\n&&\n";
        let batch = read_batch(&mut Cursor::new(r1), &mut Cursor::new(r2), 0).unwrap();
        let mut records = BamRecords::default();
        add_record(&mut records, b"one", 0x41, 1, 1);
        add_record(&mut records, b"one", 0x841, 19, 1);
        add_record(&mut records, b"one", 0x81, 1, 1);
        add_record(&mut records, b"two", 0x41, 19, 1);
        add_record(&mut records, b"two", 0x81, 1, 19);
        add_record(&mut records, b"three", 0x41, 1, 19);
        add_record(&mut records, b"three", 0x81, 1, 1);
        let (a, b, primary, other) = select(&batch, &mut records, 19).unwrap();
        assert_eq!((primary, other), (1, 2));
        assert_eq!(a, r1);
        assert_eq!(b, r2);
    }

    /// Run with the five NOMICS_SLICE_CHECK_* paths below. The FASTQs must
    /// be the 0001 fastp pair that produced the supplied golden chunk BAM.
    #[test]
    #[ignore = "requires the Spark tier1 exome chunk and GRCh38 index"]
    fn tier1_chr20_pair_names_match_golden_bam() {
        let path = |key| PathBuf::from(std::env::var(key).unwrap());
        let dir = std::env::temp_dir().join(format!("nomics-slice-check-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let out1 = dir.join("r1.fastq.gz");
        let out2 = dir.join("r2.fastq.gz");
        run(&Args {
            index: path("NOMICS_SLICE_CHECK_INDEX"),
            fastq1: path("NOMICS_SLICE_CHECK_FASTQ1"),
            fastq2: path("NOMICS_SLICE_CHECK_FASTQ2"),
            out1: out1.clone(),
            out2: out2.clone(),
            contig: "chr20".into(),
            dense_sa_cache: path("NOMICS_SLICE_CHECK_DENSE_SA"),
            dense_sa_interval: 4,
            threads: 20,
        })
        .unwrap();
        let mut bam = nomics_io::bgzf::BgzfReader::new(BufReader::new(
            File::open(path("NOMICS_SLICE_CHECK_GOLDEN_BAM")).unwrap(),
        ));
        let header = nomics_io::bam::Header::read(&mut bam).unwrap();
        let tid = header.refs.iter().position(|r| r.name == b"chr20").unwrap() as i32;
        let mut expected = BTreeSet::new();
        let mut raw = Vec::new();
        while nomics_io::bam::read_raw_record(&mut bam, &mut raw).unwrap() {
            let view = RecordView::new(&raw).unwrap();
            if view.ref_id() == tid || view.mate_ref_id() == tid {
                expected.insert(view.read_name().to_vec());
            }
        }
        for output in [&out1, &out2] {
            let mut input = Input::open(output).unwrap();
            let mut actual = BTreeSet::new();
            let mut bytes = Vec::new();
            while let Some(read) = raw_read(&mut input.reader, &mut bytes).unwrap() {
                actual.insert(read.name);
                bytes.clear();
            }
            input.finish().unwrap();
            assert_eq!(actual, expected, "{}", output.display());
        }
        std::fs::remove_dir_all(dir).unwrap();
    }
}
