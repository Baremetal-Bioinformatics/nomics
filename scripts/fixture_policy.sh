# Sourced by gates. Preserve optional fixtures by default; require them on demand.
fixture_skip() {
  echo "SKIP $*"
  if [[ ${NOMICS_REQUIRE_FIXTURES:-} == 1 ]]; then
    echo "Required fixture/oracle absent: $*" >&2
    return 1
  fi
}
