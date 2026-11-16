#!/usr/bin/env bash
# Run under lockf -k /tmp/nomics-heavy.lock; builds only this checkout's target.
# Every path uses the same two unmodified BWA -B 3 fixtures; SAM records must
# match literally, and BAMs use the existing GATES.md normalization.
# Paths: neon and portable (CPU scoring), then metal-fused (Metal seeding, fused
# seed extension and mate rescue at B=3), which fails unless its log shows the
# scoring kernels did the work and none fell back. Two-pass extension exists only
# on CUDA and is gated on the gate machine.
#   PATHS="neon portable" scripts/mac_bwa_b3_gate.sh   (CPU paths only)
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
export CARGO_TARGET_DIR="$REPO/target"
OUT=${1:-$CARGO_TARGET_DIR/b3-validation}
FIX=${NOMICS_B3_FIXTURE:-$HOME/omics-fixtures/bwa-b3-tier0}
IDX=${NOMICS_BWA_INDEX:-$HOME/omics/ref/GATK.GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64}
DENSE=${NOMICS_DENSE_SA:-$HOME/omics/cache/sa4.GATK.GRCh38.nsa}
: "${NOMICS_RG_DS:?set NOMICS_RG_DS to the read-group DS path of the goldens (see scripts/mac_tier0.sh)}"; export NOMICS_RG_DS
mkdir -p "$OUT"
for path in ${PATHS:-neon portable metal-fused}; do
    features=(); gpu_env=(NOMICS_GPU_EXTEND=fused)
    case $path in
        neon) ;;
        portable) features=(--features nomics-align/portable) ;;
        metal-fused) features=(--features nomics/metal,nomics-align/metal); gpu_env=(NOMICS_GPU_SEED=1 NOMICS_TIMING=1 NOMICS_GPU_EXTEND=fused) ;;
        *) echo "unknown path $path" >&2; exit 2 ;;
    esac
    /usr/bin/time -l cargo build --manifest-path "$REPO/Cargo.toml" --release -p nomics -p nomics-align --bin nomics --example align_chunk ${features[@]+"${features[@]}"} > "$OUT/build-$path.log" 2>&1
    for chunk in 0001 0002; do
        prefix="$OUT/$path-$chunk"
        # The CPU paths still ask for GPU extension: a build without a GPU
        # must fall back to the CPU and give the same records.
        /usr/bin/time -l env "${gpu_env[@]}" NOMICS_SAM_OUTPUT="$prefix.sam" \
            "$CARGO_TARGET_DIR/release/examples/align_chunk" "$IDX" \
            "$FIX/$chunk.HG002_WES-L1_R1.fastp.fastq.gz" \
            "$FIX/$chunk.HG002_WES-L1_R2.fastp.fastq.gz" \
            "$prefix.bam" "${NOMICS_THREADS:-8}" "$DENSE" 3 > "$prefix.log" 2>&1
        python3 - "$FIX/$chunk.B3.sam.gz" "$prefix.sam" <<'PY'
import gzip, sys
with gzip.open(sys.argv[1], 'rb') as f:
    expected = b''.join(line for line in f if not line.startswith(b'@'))
with open(sys.argv[2], 'rb') as f:
    actual = f.read()
assert actual == expected, f'SAM differs: {sys.argv[2]}'
print(f'PASS literal SAM: {sys.argv[2]}, {actual.count(bytes([10]))} records')
PY
        /usr/bin/time -l "$CARGO_TARGET_DIR/release/nomics" dev gate-bam \
            --ours "$prefix.bam" --golden "$FIX/HG002_WES.$chunk.B3.bam" > "$prefix.gate.log" 2>&1
        cat "$prefix.gate.log"
        if [[ $path == metal-* ]]; then
            if grep -E "(extension|rescue|seed search) failed|GPU FALLBACK|warm-up failed" "$prefix.log"; then
                echo "FAIL $path $chunk: a GPU stage fell back to the CPU" >&2; exit 1
            fi
            grep -qE "\[timing\] fused extend \([^)]*\): [1-9][0-9]* jobs" "$prefix.log" || { echo "FAIL $path $chunk: no GPU seed extension in the log" >&2; exit 1; }
            grep -qE "\[timing\] gpu rescue: [0-9]+ calls \([1-9][0-9]* on the GPU" "$prefix.log" || { echo "FAIL $path $chunk: no GPU mate rescue in the log" >&2; exit 1; }
            echo "PASS $path $chunk: $(grep -cE "\[timing\] fused extend" "$prefix.log") extension launches and $(grep -cE "\[timing\] gpu rescue" "$prefix.log") rescue launches on Metal, no fallback"
        fi
    done
done
