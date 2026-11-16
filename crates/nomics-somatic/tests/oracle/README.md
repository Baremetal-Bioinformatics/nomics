Numerical oracles are generated independently with OpenJDK 17 and the unmodified
Commons Math source at tag MATH_3_5 (b3c5dae8f253fcb4484e5cd3cc5662587803efc2).
They are not computed from the Rust implementation. These Java harnesses are
validation tools; neither Java nor Commons Math is a runtime dependency.

With `source` pointing at the checkout's `src/main/java` and `classes` a fresh
scratch directory outside the repository:

```sh
javac -d "$classes" -sourcepath "$source" \
  "$source/org/apache/commons/math3/util/FastMath.java" \
  "$source/org/apache/commons/math3/linear/SingularValueDecomposition.java" \
  "$source/org/apache/commons/math3/optim/univariate/BrentOptimizer.java" \
  "$source/org/apache/commons/math3/distribution/BinomialDistribution.java" \
  "$source/org/apache/commons/math3/stat/descriptive/moment/Mean.java" \
  "$source/org/apache/commons/math3/stat/descriptive/rank/Percentile.java"
javac -cp "$classes" -d "$classes" Numerics.java MoreNumerics.java
java -cp "$classes" Numerics > java_oracle_bits.txt
java -cp "$classes:$gatk_classpath" MoreNumerics > more.txt
```

`java_more_bits.txt` is the first 26 lines of `more.txt`. Its later integer and
persistence results are written explicitly in segmentation's unit tests. The
persistence class is the actual GATK 4.6.2.0 source (its only extra dependency is
`Utils` argument validation). During implementation the class was compiled in
isolation with equivalent argument checks. The reflection in MoreNumerics also
allows using it directly from a GATK runtime classpath.

The two SVD matrices exercise rank deficiency and nontrivial QR iteration; all
singular values and every U entry are checked by raw binary64 bits. Binomial
fixtures span endpoints, both deviance branches, table and asymptotic Stirling
errors, and exp underflow. Brent fixtures check point AND objective, with a
separate exhausted-budget test. Mean, LEGACY percentiles, compensated sums,
Java logarithm/rounding/formatting, Random rejection, persistence ties, and cost
recurrence also have fixtures or directly hand-evaluated arithmetic tests.

For the full-chain independent oracle, the actual GATK ContaminationModel,
ContaminationSegmenter, KernelSegmenter and PersistenceOptimizer sources were
compiled against that Commons Math build, replacing only logging, pair/interval
containers and overlap queries with equivalent validation scaffolding. The
normal and tumour filtered site counts were 147029 and 144612. The Java oracle
reproduced all 42 golden segmentation rows and the contamination/error pair.
A per-site comparison on all 144612 filtered tumour rows at c=0, maf=.4 matched
every genotype likelihood bit. It exposed 5953 differences in platform `ln`
versus Java `Math.log` (Java Math.log and StrictMath.log agreed at every site).
The Rust implementation therefore ports OpenJDK's fdlibm logarithm, rather than
using platform `ln` or substituting Commons Math FastMath.log for Java Math.log.

`HashMapOracle.java` independently checks computeIfAbsent bucket chains, String
hash collisions, red-black tree bins, root-to-front movement, and resize splitting.
Compile it with javac and run it to regenerate `hashmap_order.txt`.

The separate release CLI check, from the repository root, is:

```sh
cargo build --release -p nomics
/usr/bin/time -p bash crates/nomics-somatic/tests/oracle/cli_gate.sh
```

It leaves its fresh output directory under /tmp, recorded in
`/tmp/nomics-contamination-gate-path`, and runs all four literal cmp commands.
The integration test also covers library execution with a relocatable fixture
root and has a clean skip when that root is absent.

Full no-matched-normal execution on the published tumour pileups also ran:
contamination 0.0019024981353576948, error 2.0702953973346697E-4; the independent
Java model agreed exactly. Its segmentation matched the matched-normal run's
separate tumour model. This is additional oracle evidence, not a published
no-normal golden. Small synthetic fixtures exercise the heuristic fallback and
empty-input result separately.

Final verification on this Mac (2026-11-12):

```sh
cargo build --release -p nomics
cargo test -p nomics-somatic
cargo clippy --no-deps -p nomics-somatic -p nomics -- -D warnings
/usr/bin/time -p bash crates/nomics-somatic/tests/oracle/cli_gate.sh
NOMICS_SOMATIC_GOLDENS=/tmp/absent-nomics-somatic-goldens \
  cargo test -p nomics-somatic --test contamination_goldens contamination_goldens -- --nocapture
```

All exited zero. Seven unit tests and two integration tests passed (131.63 s for
the unoptimized integration tests). Release gate: tumour gather PASS, normal
gather PASS, contamination PASS, segmentation PASS; real 2.75 s, user 2.69 s,
system 0.03 s. All comparisons were byte-for-byte, without normalization.
Release outputs: `/tmp/nomics-contamination-gate.EUrrHz`. The absent-directory
check printed SKIP and passed. Build and Clippy emitted the existing unrelated
nomics-hc `AtomicU64::fetch_update` deprecation warning; no targeted lint failed.
Both gathered tables have zero duplicate loci.

FilterMutectCalls adds `FilteringNumerics.java` and `filtering_numerics.bits`:
70 raw-double fixtures from Commons Math 3.5 Gamma.logGamma/digamma, Beta.logBeta
and Beta.regularizedBeta. The exact JDK 17 regeneration commands are in
`docs/spec/mutect2-filtering.md` §9.3. `filtering_cli_gate.sh` runs the release CLI
with the archived inputs, literal stats cmp, complete record/header comparisons
and a BGZF-member envelope check. It never invokes a reference executable.


`FilteringLearning.java` extends the filtering oracle with 83 raw-double values:
Java scalar pow/log10/log1p/expm1, the exact SNV offset expression, an error-type
product-order witness, and three rounds of seven strand E steps separated by
two complete Commons Math Brent M steps. The fixture is
`filtering_learning.bits`; the Rust test checks every value by bits. Compile
with the same `-sourcepath` as FilteringNumerics into a fresh /tmp directory,
then run `FilteringLearning`. Full regeneration commands, all six boundary
results and the limits of the product-order witness are in
`docs/spec/mutect2-filtering.md` §9.4–§9.5. The five-case integration root is
`NOMICS_SOMATIC_FMC_PERTURB`, default `~/omics/runs/somatic/fmc-perturb`; missing
fixtures print SKIP. Neither Java harness is used at runtime.
