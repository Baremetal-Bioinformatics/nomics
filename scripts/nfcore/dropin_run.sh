#!/usr/bin/env bash
# nf-core/sarek under Nextflow with the nomics shims on PATH in every task ("drop-in"), then
# gated against a stock sarek run of the same input (docs/strategy/nfcore-dropin.md, phase 2).
# On the Spark. reference/run_sarek.sh takes the resource lock itself: do not wrap this in flock.
#
#   scripts/nfcore/dropin_run.sh tier0|tier1 [RUN_DIR]
#
# 1. builds a CPU-only nomics (no CUDA in the task containers) and the shim directory
# 2. runs sarek with reference/nomics-dropin.config
# 3. shim gate: every task of a process the shims are meant to handle logged
#    `nomics-shim: handled`; a delegated or silent one fails
# 4. output gate: reference/compare_runs.sh against the golden run (records, headers, tables),
#    and the reports/ tree (bytes; fastp under its exempt-field rule; metrics but lines 2 and 4)
set -uo pipefail
tier=${1:?usage: dropin_run.sh tier0|tier1 [RUN_DIR]}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
case $tier in
  tier0) SHEET=$REPO/reference/samplesheet.tier0.csv; GOLD=$HOME/omics/runs/tier0-wes-20261008-r1 ;;
  tier1) SHEET=$REPO/reference/samplesheet.tier1.csv; GOLD=$HOME/omics/runs/tier1-wes-20261008-r1 ;;
  *) echo "tier must be tier0 or tier1"; exit 2 ;;
esac
RUN=${2:-$HOME/omics/runs/dropin-$tier-$(date +%Y%m%d-%H%M)}
export PATH=$HOME/.cargo/bin:$PATH
export NOMICS_SHIM_BIN=${NOMICS_SHIM_BIN:-$HOME/omics/nomics-shim/bin}

CARGO_TARGET_DIR=$HOME/omics/build/cpu-dropin cargo build --manifest-path "$REPO/Cargo.toml" --release -p nomics 2>&1 | tail -n 1
mkdir -p "$NOMICS_SHIM_BIN"
cp "$HOME/omics/build/cpu-dropin/release/nomics" "$NOMICS_SHIM_BIN/nomics.new" && mv "$NOMICS_SHIM_BIN/nomics.new" "$NOMICS_SHIM_BIN/nomics"
for t in fastp bwa gatk samtools; do ln -sfn nomics "$NOMICS_SHIM_BIN/$t"; done

# DROPIN_GATE_ONLY=1 gates an existing RUN_DIR again without rerunning sarek.
[[ -n ${DROPIN_GATE_ONLY:-} ]] || NOMICS_SAREK_CONFIG=$REPO/reference/nomics-dropin.config "$REPO/reference/run_sarek.sh" "$RUN" "$SHEET" "$REPO/reference/params.tier1.yaml" > "$RUN.log" 2>&1
echo "sarek: $(cat "$RUN/exit_status" 2>/dev/null), $(grep -m1 'Elapsed (wall clock)' "$RUN/time.txt" | sed 's/^[[:space:]]*//')"
[[ $(cat "$RUN/exit_status" 2>/dev/null) == EXIT=0 ]] || { tail -n 30 "$RUN.log"; exit 1; }

# ---- shim gate: what ran as nomics, per process
fail=0
python3 -I - "$RUN" <<'PY' || fail=1
import collections, csv, os, re, sys
run = sys.argv[1]
expected = {'FASTP': 'fastp', 'BWAMEM1_MEM': 'bwa mem', 'GATK4_MARKDUPLICATES': 'gatk MarkDuplicates',
            'GATK4_BASERECALIBRATOR': 'gatk BaseRecalibrator', 'GATK4_APPLYBQSR': 'gatk ApplyBQSR',
            'GATK4_HAPLOTYPECALLER': 'gatk HaplotypeCaller'}
info = os.path.join(run, 'results', 'pipeline_info')
trace = sorted(f for f in os.listdir(info) if f.startswith('execution_trace'))
rows = list(csv.DictReader(open(os.path.join(info, trace[-1])), delimiter='\t'))
handled, delegated, silent, bad = collections.Counter(), collections.Counter(), collections.Counter(), []
for r in rows:
    proc = r['name'].split(' ')[0].split(':')[-1]
    h = r['hash'].split('/')
    d = [x for x in os.listdir(os.path.join(run, 'work', h[0])) if x.startswith(h[1])]
    err = open(os.path.join(run, 'work', h[0], d[0], '.command.err'), errors='replace').read() if d else ''
    lines = re.findall(r'nomics-shim: ((?:handled|delegated) [^\n]*?)(?=nomics-shim: |$)', err, flags=re.M)
    for l in lines:
        (handled if l.startswith('handled') else delegated)[(proc, l)] += 1
    if not lines:
        silent[proc] += 1
    if proc in expected and not any(l == 'handled ' + expected[proc] for l in lines):
        bad.append((proc, r['hash'], lines))
print('handled:');   [print(f'  {n:3d}  {p}: {l}') for (p, l), n in sorted(handled.items())]
print('delegated:'); [print(f'  {n:3d}  {p}: {l}') for (p, l), n in sorted(delegated.items())]
print('no shim call:', ', '.join(f'{p} x{n}' for p, n in sorted(silent.items())))
for b in bad: print('NOT HANDLED:', b)
print(f'shim gate: {sum(1 for r in rows if r["name"].split(" ")[0].split(":")[-1] in expected) - len(bad)} expected tasks handled, {len(bad)} not')
sys.exit(1 if bad else 0)
PY

# ---- output gate
"$REPO/reference/compare_runs.sh" "$GOLD" "$RUN" > "$RUN/compare.tsv" 2> "$RUN/compare.err" || { echo "compare_runs failed: $(tail -n 2 "$RUN/compare.err")"; fail=1; }
# A file passes when its strictest passing check is one of: bytes, records+header (BAM, VCF), body (table).
python3 -I - "$RUN/compare.tsv" <<'PY' || fail=1
import collections, sys
checks = collections.defaultdict(dict)
for line in open(sys.argv[1]):
    stage, f, check, result, detail = (line.rstrip('\n').split('\t') + [''] * 5)[:5]
    checks[f][check] = (result, detail)
bad = 0
for f, c in sorted(checks.items()):
    if c.get('bytes', ('',))[0] == 'PASS': verdict = 'bytes'
    elif f.endswith('.sorted.bam'): verdict = 'not gated: samtools merge gives its read groups random suffixes in every run (docs/GATES.md)'
    elif c.get('records', ('',))[0] == 'PASS' and c.get('header', ('',))[0] == 'PASS': verdict = 'records + header'
    elif c.get('body', ('',))[0] == 'PASS': verdict = 'table body'
    else: verdict = None
    if verdict is None:
        bad += 1; print('FAIL', f, {k: v for k, v in c.items() if v[0] != 'PASS'})
    else: print('PASS', f, f'({verdict})')
print(f'output gate: {len(checks) - bad} pass, {bad} fail')
sys.exit(1 if bad else 0)
PY

norm() { python3 -I "$REPO/scripts/fastp_norm.py" "$1" "$2"; }
rp=0 rf=0
while IFS= read -r g; do
  rel=${g#"$GOLD"/results/reports/}; o=$RUN/results/reports/$rel
  [[ -f $o ]] || { rf=$((rf + 1)); echo "FAIL reports/$rel (missing)"; continue; }
  case $rel in
    fastp/*.fastp.json|fastp/*.fastp.html|fastp/*.fastp.log) e=${rel##*.}; cmp -s <(norm "$e" "$o") <(norm "$e" "$g") ;;
    markduplicates/*.metrics) cmp -s <(sed '2d;4d' "$o") <(sed '2d;4d' "$g") ;;
    fastqc/*.zip) cmp -s <(unzip -p "$o" '*/fastqc_data.txt' '*/summary.txt') <(unzip -p "$g" '*/fastqc_data.txt' '*/summary.txt') ;;
    fastqc/*.html) true ;;   # rendered by the real FastQC in both runs; carries no gate
    mosdepth/*.gz|mosdepth/*.csi) cmp -s <(gzip -dc "$o") <(gzip -dc "$g") ;;
    vcftools/*.TsTv.qual)
      # vcftools leaves stray values in some fields, different in every reference run (docs/GATES.md):
      # identical, or both sides pass nomics' exempt-field comparison against the run's own VCF.
      cmp -s "$o" "$g" || {
        vt=$(dirname "$o"); sub=${vt#"$RUN"/results/reports/vcftools/}; base=$(basename "$o" .TsTv.qual)
        [[ $("$NOMICS_SHIM_BIN/nomics" qc vcf-gate --vcf "$RUN/results/variant_calling/$sub/$base.vcf.gz" --vcftools "$vt" --bcftools "$RUN/results/reports/bcftools/$sub" 2>&1) == *"verdict: PASS"* ]]
      } ;;
    *) cmp -s "$o" "$g" ;;
  esac && rp=$((rp + 1)) || { rf=$((rf + 1)); echo "FAIL reports/$rel"; }
done < <(find -L "$GOLD/results/reports" -type f | sort)
echo "reports gate: $rp pass, $rf fail"
(( rf == 0 )) || fail=1
echo "DROPIN_DONE $tier fail=$fail  ($RUN)"
exit $fail
