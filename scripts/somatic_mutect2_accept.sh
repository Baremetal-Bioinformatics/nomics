#!/usr/bin/env bash
# Lead acceptance of the Mutect2 caller on the Spark: every caller gate on the ten slices of
# the HCC1395 exome pair, repeat runs with peak memory, then the whole pair through the chain
# gate (scripts/somatic_chain_gate.sh). About 100 minutes; run it under the resource lock.
#
#   scripts/somatic_mutect2_accept.sh [OUT_DIR]      NOMICS=path/to/nomics
#   defaults: OUT_DIR=~/omics/runs/somatic/accept-<short SHA or "tree">, NOMICS=target/release/nomics
#
# The slices are ~/omics/runs/somatic/m2-<name> (scripts/somatic_mutect2_slice.sh); chr17 and
# chr11 were held out from the port. The internal haplotype and read-tag checks on chr21 need
# the Java oracle dumps (NOMICS_M2_ASSEMBLY_ORACLE, NOMICS_M2_LIKELIHOODS_ORACLE); without
# them those two checks print SKIP, and this script counts and reports the skips.
# Linux (GNU time, docker for the GATK image's java). Rules: docs/GATES.md.
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd); cd "$REPO" || exit 2
export NOMICS=${NOMICS:-${CARGO_TARGET_DIR:-$REPO/target}/release/nomics}
[[ -x $NOMICS ]] || { echo "no nomics binary at $NOMICS (cargo build --release -p nomics)"; exit 2; }
S=$HOME/omics/runs/somatic
OUT=${1:-$S/accept-$(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo tree)}
rm -rf "$OUT"; mkdir -p "$OUT/bin"
REF=${REF:-$HOME/omics/ref/GATK.GRCh38}; FA=$REF/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
IMG=${GATK_IMAGE:-community.wave.seqera.io/library/gatk4_htslib_samtools:6fd46722be032bdd}
for t in java javac; do
  printf '#!/usr/bin/env bash\nexec docker run --rm -i -u %s:%s -v %s/omics:%s/omics -w "$PWD" %s %s "$@"\n' \
    "$(id -u)" "$(id -g)" "$HOME" "$HOME" "$IMG" "$t" > "$OUT/bin/$t"; chmod +x "$OUT/bin/$t"
done
export JAVA=$OUT/bin/java JAVAC=$OUT/bin/javac
ALL=${ACCEPT_SLICES:-"chr17 chr11 chr21 chr22 chr19 chr20 chr13 chr16 chr12 multi"}
fails=0; skips=0
step() { # step NAME EXIT_STATUS
  echo "exit=$2"; (( $2 == 0 )) || { fails=$((fails+1)); echo "FAILED: $1"; }
}
f() { grep -E "records identical|header:|stats:|f1r2|priors|^real|FAIL|SKIP|\| [0-9]+ \|"; }

echo "== nomics $("$NOMICS" --version 2>/dev/null | head -n 1); tree $(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo "not a git checkout"); $(date -u +%FT%TZ)"
echo "== final gate, ten slices"
NOMICS_M2_CALLS_OUT=$OUT/calls bash scripts/somatic_mutect2_gate.sh $ALL | tee "$OUT/final.log" | f; step "final gate" "${PIPESTATUS[0]}"
echo "== regions gate"
bash scripts/somatic_mutect2_regions_gate.sh $ALL | tee "$OUT/regions.log"; step "regions gate" "${PIPESTATUS[0]}"
echo "== standalone f1r2 gate"
M2_F1R2_OUT=$OUT/f1r2 bash scripts/somatic_mutect2_f1r2_gate.sh $ALL | tee "$OUT/f1r2.log"; step "f1r2 gate" "${PIPESTATUS[0]}"
echo "== haplotypes gate"
bash scripts/somatic_mutect2_haplotypes_gate.sh "$OUT/hap" $ALL > "$OUT/hap.log" 2>&1; rc=$?
grep -E "PASS|FAIL|SKIP" "$OUT/hap.log"; step "haplotypes gate" "$rc"
echo "== read-tag gate"
bash scripts/somatic_mutect2_likelihoods_gate.sh "$OUT/xa" $ALL > "$OUT/xa.log" 2>&1; rc=$?
grep -E "PASS|FAIL|SKIP" "$OUT/xa.log"; step "read-tag gate" "$rc"
rm -rf "$OUT/xa" "$OUT/hap"
skips=$(cat "$OUT"/*.log | grep -c "^SKIP" || true)

echo "== repeats and peak memory"
for c in chr17 chr19 multi; do d=$S/m2-$c
  [[ " $ALL " == *" $c "* ]] || continue
  /usr/bin/time -v "$NOMICS" somatic mutect2 --input "$d/HCC1395BL_WES.$c.bam" --input "$d/HCC1395_WES.$c.bam" \
    --normal-sample HCC1395_HCC1395BL_WES --reference "$FA" --intervals "$d/$c.bed" \
    --panel-of-normals "$d/1000g_pon.$c.vcf.gz" --germline-resource "$d/af-only-gnomad.$c.vcf.gz" \
    --output "$OUT/rep.$c.vcf.gz" --f1r2-tar-gz "$OUT/rep.$c.f1r2.tar.gz" 2>&1 | grep -E "Maximum resident|Elapsed" | sed "s/^/$c /"
  if cmp <(gzip -dc "$OUT/rep.$c.vcf.gz") <(gzip -dc "$OUT/calls/$c.vcf.gz") && cmp "$OUT/rep.$c.f1r2.tar.gz" "$OUT/calls/$c.f1r2.tar.gz"
  then echo "$c repeat: VCF and archive identical"; else fails=$((fails+1)); echo "FAILED: $c repeat differs"; fi
done

if [[ -z ${ACCEPT_NO_CHAIN:-} ]]; then
  echo "== whole pair: chain gate $(date -u +%FT%TZ)"
  bash scripts/somatic_chain_gate.sh "$S/S1-r1" "$OUT/chain"; step "chain gate" "$?"
fi
echo "== done $(date -u +%FT%TZ)"
echo "somatic acceptance: $fails failed step(s), $skips skipped check(s) (outputs in $OUT)"
(( fails == 0 ))
