#!/usr/bin/env python3
"""Independent BGZF check with Python's zlib (not nomics' inflate).

    scripts/bgzf_check.py file.bam [more.bam ...]

For every member: the gzip/BGZF header (magic, FEXTRA, the BC subfield,
BSIZE), raw-DEFLATE decompression with zlib (wbits=-15) consuming exactly the
compressed data, then CRC32 and ISIZE against the footer. The file must end
with the 28-byte EOF member. Prints block counts, payload-size histogram of
non-full blocks, and the SHA-256 of the uncompressed stream; exits 1 on any
error.
"""
import hashlib
import struct
import sys
import zlib

EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
FULL = 65498


def check(path):
    sha = hashlib.sha256()
    blocks = full = 0
    partial = {}
    total_u = 0
    with open(path, "rb") as f:
        data_at = 0
        last = b""
        while True:
            head = f.read(18)
            if not head:
                break
            if len(head) < 18:
                raise ValueError(f"{path}: truncated header at {data_at}")
            id1, id2, cm, flg, _mtime, _xfl, _os, xlen = struct.unpack("<BBBBIBBH", head[:12])
            si1, si2, slen, bsize = struct.unpack("<BBHH", head[12:18])
            if (id1, id2, cm) != (0x1F, 0x8B, 8) or not flg & 4 or xlen != 6 or (si1, si2, slen) != (66, 67, 2):
                raise ValueError(f"{path}: bad BGZF header at {data_at}")
            rest = f.read(bsize + 1 - 18)
            if len(rest) != bsize + 1 - 18:
                raise ValueError(f"{path}: truncated member at {data_at}")
            cdata, crc, isize = rest[:-8], *struct.unpack("<II", rest[-8:])
            d = zlib.decompressobj(-15)
            payload = d.decompress(cdata)
            if not d.eof or d.unused_data or d.unconsumed_tail:
                raise ValueError(f"{path}: deflate stream does not end exactly at member end ({data_at})")
            if len(payload) != isize or zlib.crc32(payload) != crc:
                raise ValueError(f"{path}: CRC32/ISIZE mismatch at {data_at}")
            sha.update(payload)
            blocks += 1
            total_u += isize
            if isize == FULL:
                full += 1
            else:
                partial[isize] = partial.get(isize, 0) + 1
            last = head + rest
            data_at += bsize + 1
    if last != EOF:
        raise ValueError(f"{path}: does not end with the BGZF EOF member")
    print(f"{path}: OK {blocks} members ({full} full of {FULL} B), {total_u} B uncompressed, "
          f"compressed {data_at} B; other payload sizes {dict(sorted(partial.items()))}; sha256 {sha.hexdigest()}")


def main():
    bad = 0
    for p in sys.argv[1:]:
        try:
            check(p)
        except ValueError as e:
            print(f"FAIL {e}")
            bad += 1
    sys.exit(1 if bad or len(sys.argv) < 2 else 0)


if __name__ == "__main__":
    main()
