//! Exact downsampler request tape. No read payload survives the plan pass.
use super::*;
use std::io::BufReader;

pub(super) struct Generator {
    pub rng: JavaRandom,
    pub calls: u64,
    tape: Option<BufWriter<File>>,
}
impl Generator {
    pub fn new(rng: JavaRandom) -> Self {
        Self {
            rng,
            calls: 0,
            tape: None,
        }
    }
    pub fn next_int(&mut self, bound: i32) -> io::Result<i32> {
        if let Some(tape) = &mut self.tape {
            tape.write_all(&bound.to_le_bytes())?;
        }
        self.calls += 1;
        Ok(self.rng.next_int(bound))
    }
}
/// Result of the front traversal, before the group-ordered replay.
pub struct Tape {
    pub contig: usize,
    pub calls: u64,
    ordinal: u32,
    path: PathBuf,
}
/// Exact whole-generator snapshots, including the Gaussian cache.
#[derive(Clone, Debug)]
pub struct Start {
    pub contig: usize,
    pub calls: u64,
    pub(crate) rng: JavaRandom,
    pub(crate) ordinal: u32,
    pub(crate) end_rng: JavaRandom,
    pub(crate) end_ordinal: u32,
}
/// Execute the production loader, post-transform, pool and reservoir path.
/// Reservoir choices cannot alter subsequent requests: the loop walks every
/// MQ>50 input in each pool irrespective of which reservoir slot was chosen.
pub fn record(inputs: &Inputs<'_>, contig: usize, path: &Path) -> io::Result<Tape> {
    let header = BamQuery::open(&inputs.input[0], &[])?.header;
    let names: Vec<_> = header
        .refs
        .iter()
        .map(|r| (r.name.clone(), i64::from(r.len)))
        .collect();
    let mut local = Intervals::from_bed(inputs.intervals, &names)?;
    if contig >= names.len() || local.unpadded[contig].is_empty() {
        return Err(io::Error::other("plan requires an occupied contig"));
    }
    for c in 0..names.len() {
        if c != contig {
            local.padded[c].clear();
            local.unpadded[c].clear();
        }
    }
    let reference = reference(inputs.reference, &names, &local)?;
    let mut stream = ReadStream {
        loader: reads::Loader::new(inputs.input, &local, inputs.normal_sample, 0)?,
        intervals: &local,
        reference: &reference,
        normal: inputs.normal_sample,
        head: None,
        pool: Vec::new(),
        stats: Stats::default(),
    };
    let mut rng = Generator::new(JavaRandom::gatk());
    rng.tape = Some(BufWriter::new(File::create(path)?));
    let (mut all, mut ids) = (Vec::new(), Vec::new());
    while stream.append(&mut all, &mut ids, &mut rng)? {
        all.clear();
        ids.clear();
    }
    rng.tape.as_mut().unwrap().flush()?;
    Ok(Tape {
        contig,
        calls: rng.calls,
        ordinal: stream.loader.ordinal(),
        path: path.into(),
    })
}
/// Replay actual nextInt calls, including rejection sampling, in dictionary order.
pub fn replay(tapes: &[Tape]) -> io::Result<Vec<Start>> {
    let mut rng = JavaRandom::gatk();
    let mut ordinal = 0u32;
    let mut starts = Vec::new();
    for tape in tapes {
        if starts
            .last()
            .is_some_and(|p: &Start| p.contig >= tape.contig)
        {
            return Err(io::Error::other("plan tapes must be in dictionary order"));
        }
        let start_rng = rng.clone();
        let mut reader = BufReader::new(File::open(&tape.path)?);
        if reader.get_ref().metadata()?.len() != tape.calls * 4 {
            return Err(io::Error::other("truncated random request tape"));
        }
        for _ in 0..tape.calls {
            let mut bytes = [0; 4];
            reader.read_exact(&mut bytes)?;
            let bound = i32::from_le_bytes(bytes);
            if bound <= 0 {
                return Err(io::Error::other("invalid random request"));
            }
            rng.next_int(bound);
        }
        let end_ordinal = ordinal.wrapping_add(tape.ordinal);
        starts.push(Start {
            contig: tape.contig,
            calls: tape.calls,
            rng: start_rng,
            ordinal,
            end_rng: rng.clone(),
            end_ordinal,
        });
        ordinal = end_ordinal;
    }
    Ok(starts)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn replay_rejection_sampling_and_wrapping_ordinals() {
        let path = std::env::temp_dir().join(format!("nomics-rng-tape-{}", std::process::id()));
        let bounds = [51i32, 64, 1_073_741_825, 2_000_000_000, 53];
        std::fs::write(
            &path,
            bounds
                .iter()
                .flat_map(|n| n.to_le_bytes())
                .collect::<Vec<_>>(),
        )
        .unwrap();
        let tapes = [
            Tape {
                contig: 1,
                calls: 5,
                ordinal: u32::MAX,
                path: path.clone(),
            },
            Tape {
                contig: 4,
                calls: 5,
                ordinal: 7,
                path: path.clone(),
            },
        ];
        let starts = replay(&tapes).unwrap();
        let mut serial = JavaRandom::gatk();
        for s in &starts {
            assert_eq!(serial, s.rng);
            for b in bounds {
                serial.next_int(b);
            }
            assert_eq!(serial, s.end_rng);
        }
        assert_eq!(starts[1].ordinal, u32::MAX);
        assert_eq!(starts[1].end_ordinal, 6);
        std::fs::remove_file(path).unwrap();
    }
}
