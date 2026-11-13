"""Generate zlib reference outputs for nomics-io deflate tests.

Run on a host whose zlib matches the reference (classic zlib 1.3.x, e.g. the Spark):
    python3 -I gen_zlib_fixtures.py <out_dir>
Inputs come from `gen()`, mirrored exactly by `synth` in tests/deflate_zlib.rs.
"""
import ctypes, ctypes.util, os, sys, zlib

MASK = (1 << 64) - 1

def gen(seed, n, alphabet, copies=True):
    """xorshift64* stream; bytes from `alphabet`, with back-copies (if `copies`) to create matches."""
    x = (seed * 0x9E3779B97F4A7C15 + 1) & MASK
    def nxt():
        nonlocal x
        x ^= x >> 12; x ^= (x << 25) & MASK; x ^= x >> 27
        return (x * 0x2545F4914F6CDD1D) & MASK
    out = bytearray()
    while len(out) < n:
        r = nxt()
        if copies and len(out) > 16 and (r & 7) == 0:
            dist = 1 + (nxt() % min(len(out), 40000))
            length = 3 + (nxt() % 300)
            for _ in range(length):
                if len(out) >= n: break
                out.append(out[-dist])
        else:
            out.append(alphabet[(r >> 8) % len(alphabet)])
    return bytes(out)

DNA = b"ACGT"
QUAL = bytes(range(33, 75))
ALL = bytes(range(256))

CASES = {
    "empty": b"",
    "one": b"x",
    "hello": b"hello",
    "dna_1k": gen(1, 1000, DNA),
    "qual_20k": gen(2, 20000, QUAL),
    "random_5k": gen(3, 5000, ALL),
    "random_65498": gen(4, 65498, ALL),
    "dna_65498": gen(5, 65498, DNA + b"N"),
    "mixed_65498": gen(6, 30000, QUAL) + gen(7, 35498, DNA),
    "zeros_65498": bytes(65498),
    "noise_65498": gen(8, 65498, ALL, copies=False),
    "noise_40000": gen(9, 40000, ALL, copies=False),
}

def raw(data, level):
    c = zlib.compressobj(level, zlib.DEFLATED, -15, 8, zlib.Z_DEFAULT_STRATEGY)
    return c.compress(data) + c.flush()

class ZStream(ctypes.Structure):
    _fields_ = [("next_in", ctypes.c_void_p), ("avail_in", ctypes.c_uint), ("total_in", ctypes.c_ulong),
                ("next_out", ctypes.c_void_p), ("avail_out", ctypes.c_uint), ("total_out", ctypes.c_ulong),
                ("msg", ctypes.c_char_p), ("state", ctypes.c_void_p), ("zalloc", ctypes.c_void_p),
                ("zfree", ctypes.c_void_p), ("opaque", ctypes.c_void_p), ("data_type", ctypes.c_int),
                ("adler", ctypes.c_ulong), ("reserved", ctypes.c_ulong)]

def htsjdk_sequence(blocks, level=2, capacity=65518):
    """One deflater, deflateReset per block, one deflate(Z_FINISH) into `capacity` bytes."""
    libz = ctypes.CDLL(ctypes.util.find_library("z"))
    libz.zlibVersion.restype = ctypes.c_char_p
    s = ZStream()
    ver = libz.zlibVersion()
    assert libz.deflateInit2_(ctypes.byref(s), level, 8, -15, 8, 0, ver, ctypes.sizeof(ZStream)) == 0
    outs = []
    for b in blocks:
        assert libz.deflateReset(ctypes.byref(s)) == 0
        inbuf = ctypes.create_string_buffer(b, len(b))
        outbuf = ctypes.create_string_buffer(capacity)
        s.next_in, s.avail_in = ctypes.addressof(inbuf), len(b)
        s.next_out, s.avail_out = ctypes.addressof(outbuf), capacity
        rc = libz.deflate(ctypes.byref(s), 4)  # Z_FINISH
        n = capacity - s.avail_out
        outs.append((rc == 1, outbuf.raw[:n]))  # rc 1 = Z_STREAM_END
    libz.deflateEnd(ctypes.byref(s))
    return outs, ver.decode()

def main(out):
    os.makedirs(out, exist_ok=True)
    for name, data in CASES.items():
        for level in (1, 2, 3):
            with open(os.path.join(out, f"{name}.l{level}.deflate"), "wb") as f:
                f.write(raw(data, level))
    # Stateful sequence: later blocks see stale window bytes from earlier ones.
    seq = [CASES["dna_65498"], CASES["mixed_65498"], CASES["dna_1k"], CASES["noise_65498"],
           CASES["random_65498"], CASES["qual_20k"], CASES["noise_40000"], CASES["dna_65498"]]
    outs, ver = htsjdk_sequence(seq)
    for i, (finished, data) in enumerate(outs):
        with open(os.path.join(out, f"seq{i}.{'done' if finished else 'outofspace'}.deflate"), "wb") as f:
            f.write(data)
    with open(os.path.join(out, "ZLIB_VERSION"), "w") as f:
        f.write(f"{ver} (python zlib runtime {zlib.ZLIB_RUNTIME_VERSION})\n")
    print("zlib", ver)

def main_lazy(out):
    """Levels 4-9 (deflate_slow) only, so the level 1-3 and sequence fixtures stay as made."""
    os.makedirs(out, exist_ok=True)
    for name, data in CASES.items():
        for level in (4, 5, 6, 9):
            with open(os.path.join(out, f"{name}.l{level}.deflate"), "wb") as f:
                f.write(raw(data, level))
    with open(os.path.join(out, "ZLIB_VERSION_LAZY"), "w") as f:
        f.write(f"levels 4-9: python zlib {zlib.ZLIB_VERSION} (runtime {zlib.ZLIB_RUNTIME_VERSION})\n")

if __name__ == "__main__":
    # --lazy: deflate_slow fixtures (levels 4, 5, 6, 9) only.
    if sys.argv[1] == "--lazy":
        main_lazy(sys.argv[2])
    else:
        main(sys.argv[1])
