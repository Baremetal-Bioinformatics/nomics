#!/usr/bin/env python3
"""Print a fastp report with the fields docs/GATES.md exempts blanked, and nothing else.

    fastp_norm.py json|html|log FILE

json: the number at duplication.rate. html: the two timestamps and the duplication rate.
log: the duplication rate and the elapsed seconds. A field that is not found is an error.
"""
import re
import sys

kind, path = sys.argv[1], sys.argv[2]
s = open(path, 'rb').read()


def one(pattern, text):
    out, n = re.subn(pattern, rb'\g<1>X\g<2>', text, count=1, flags=re.S)
    if n != 1:
        sys.exit(f'{path}: field not found: {pattern!r}')
    return out


if kind == 'json':
    s = one(rb'("duplication": \{\s*"rate": )[-+.0-9eEna]+(\s*\})', s)
elif kind == 'html':
    s = one(rb'(<title>fastp report at )[^<]*?( </title>)', s)
    s = one(rb"(duplication rate:</td><td class='col2'>)[-+.0-9eEna]+(%</td>)", s)
    s = one(rb'(fastp [0-9.]+, at )[^<]*?( </div></body></html>)\Z', s)
elif kind == 'log':
    s = one(rb'(\nDuplication rate: )[-+.0-9eEna]+(%\n)', s)
    s = one(rb'(, time used: )[0-9]+( seconds\n)', s)
else:
    sys.exit(f'unknown kind {kind!r}')
sys.stdout.buffer.write(s)
