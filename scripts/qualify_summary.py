#!/usr/bin/env python3
"""Summarize the tab-separated NOMICS_FIXTURE_LOG (stdlib only)."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path


def summarize(path):
    counts = defaultdict(Counter)
    skipped = []
    data = Path(path).read_text()
    if data and not data.endswith('\n'):
        raise ValueError(f'{path}: unterminated fixture log record')
    for number, line in enumerate(data.splitlines(), 1):
        fields = line.split('\t')
        if len(fields) != 4 or fields[2] not in ('executed', 'skipped'):
            raise ValueError(f'{path}:{number}: invalid fixture log record')
        test, crate, status, wanted = fields
        counts[crate][status] += 1
        if status == 'skipped':
            skipped.append((crate, test, wanted))
    return counts, skipped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('log')
    parser.add_argument('--require-all', action='store_true')
    args = parser.parse_args()
    try:
        counts, skipped = summarize(args.log)
    except (OSError, ValueError) as error:
        parser.exit(2, f'{error}\n')
    total = Counter()
    for crate, count in sorted(counts.items()):
        print(f'{crate}: executed={count["executed"]} skipped={count["skipped"]}')
        total.update(count)
    print(f'TOTAL: executed={total["executed"]} skipped={total["skipped"]}')
    for crate, test, wanted in skipped:
        print(f'SKIPPED {crate}::{test}: {wanted}')
    return int(total['executed'] == 0 or (args.require_all and bool(skipped)))


if __name__ == '__main__':
    raise SystemExit(main())
