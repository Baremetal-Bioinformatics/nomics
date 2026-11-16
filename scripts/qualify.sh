#!/usr/bin/env bash
# Run on the gate machine with its fixture/oracle environment configured.
set -uo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo" || exit 1
# Where this machine keeps its fixtures and oracles: a file of VAR=value lines kept
# outside the repository (template: scripts/qualify.env.example). Read before anything
# else so it can set TMPDIR and the test profile; it cannot turn the requirement off.
envfile=${NOMICS_QUALIFY_ENV:-}
if [[ -n $envfile ]]; then
  [[ -r $envfile ]] || { echo "error: NOMICS_QUALIFY_ENV is set but not readable: $envfile" >&2; exit 1; }
  set -a
  # shellcheck disable=SC1090
  . "$envfile" || exit 1
  set +a
fi
out=${1:-$(mktemp -d "${TMPDIR:-/tmp}/nomics-qualify.XXXXXX")}
mkdir -p "$out" || exit 1
out=$(cd "$out" && pwd)
export CARGO_TARGET_DIR="$repo/target"
export NOMICS_REQUIRE_FIXTURES=1 NOMICS_FIXTURE_LOG="$out/fixtures.tsv"
: > "$NOMICS_FIXTURE_LOG"
# Override a caller's CLI selection: qualify the release binary built here.
export NOMICS_QC_BIN="$repo/target/release/nomics"
set --
[[ -z ${NOMICS_QUALIFY_FEATURES:-} ]] || set -- --features "$NOMICS_QUALIFY_FEATURES"
build=0 tests=0 summary=0
cargo build --release -p nomics "$@" > "$out/build.log" 2>&1 || build=$?
# Keep test failures and missing-fixture evidence even when the release build fails.
cargo test --workspace --no-fail-fast "$@" > "$out/tests.log" 2>&1 || tests=$?
python3 -I scripts/qualify_summary.py "$NOMICS_FIXTURE_LOG" --require-all > "$out/summary.txt" 2>&1 || summary=$?
{
  # A tree copied to the gate machine without its .git names its commit in NOMICS_TREE_SHA.
  if git rev-parse HEAD > /dev/null 2>&1; then
    echo "git: $(git rev-parse HEAD)"
    echo "working tree: $(git status --porcelain | wc -l | tr -d ' ') changed paths"
  else
    echo "git: ${NOMICS_TREE_SHA:-unknown (not a git checkout; set NOMICS_TREE_SHA)}"
  fi
  rustc --version
  cargo --version
  echo "host: $(hostname)"
  if [[ -r /proc/cpuinfo ]]; then
    cpu=$(awk -F ': ' '/model name/ { print $2; exit }' /proc/cpuinfo)
    if [[ -z $cpu ]] && command -v lscpu >/dev/null; then
      cpu=$(LC_ALL=C lscpu | sed -n 's/^Model name: *//p')
    fi
    echo "CPU: ${cpu:-unknown (see /proc/cpuinfo)}"
  else
    echo "CPU: $(sysctl -n machdep.cpu.brand_string)"
  fi
  echo "features: ${NOMICS_QUALIFY_FEATURES:-default}"
  echo "test profile: opt-level ${CARGO_PROFILE_DEV_OPT_LEVEL:-default}, overflow checks ${CARGO_PROFILE_DEV_OVERFLOW_CHECKS:-default}, debug assertions ${CARGO_PROFILE_DEV_DEBUG_ASSERTIONS:-default}"
  if [[ -n $envfile ]]; then
    if command -v sha256sum >/dev/null; then digest=$(sha256sum < "$envfile"); else digest=$(shasum -a 256 < "$envfile"); fi
    echo "environment file: $(basename "$envfile"), sha256 ${digest%% *}, $(grep -cE '^[[:space:]]*(export[[:space:]]+)?[A-Za-z_][A-Za-z0-9_]*=' "$envfile") variables"
  else
    echo "environment file: none"
  fi
  echo "exit codes: build=$build tests=$tests summary=$summary"
  cat "$out/summary.txt"
} > "$out/report.txt"
cat "$out/report.txt"
echo "Evidence: $out"
(( build == 0 && tests == 0 && summary == 0 ))
