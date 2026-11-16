// Oracle for nomics-hc's PairHMM: GATK 4.6.2.0 LoglessPairHMM on generated cases.
//   java -cp gatk-package-4.6.2.0-local.jar GatkPairHmm.java N SEED > pairhmm.tsv
// Lines: "hap read quals ins del gcp hexResult" (quals as comma-separated ints).
// Cases mimic HC: haplotypes are mutated random sequences, reads are noisy substrings,
// base quals are 6 or 18..41 (after capMinimumReadQualities), ins/del quals 6..45, gcp 10;
// one HMM instance is reused across a batch (as GATK does per region).
package org.broadinstitute.hellbender.utils.pairhmm;

import java.util.*;

public class GatkPairHmm {
    static final byte[] ACGT = {'A', 'C', 'G', 'T'};

    static String join(byte[] q) {
        StringBuilder sb = new StringBuilder();
        for (byte b : q) { if (sb.length() > 0) sb.append(','); sb.append(b); }
        return sb.toString();
    }

    public static void main(String[] a) {
        int n = Integer.parseInt(a[0]);
        Random r = new Random(Long.parseLong(a[1]));
        StringBuilder out = new StringBuilder();
        int done = 0;
        while (done < n) {
            // One "region": a few haplotypes and reads, HMM initialized to the batch maxima.
            int hapLen = 20 + r.nextInt(430);
            byte[] base = new byte[hapLen];
            for (int i = 0; i < hapLen; ++i) base[i] = ACGT[r.nextInt(4)];
            List<byte[]> haps = new ArrayList<>();
            int nh = 1 + r.nextInt(4);
            for (int h = 0; h < nh; ++h) {
                byte[] hb = base.clone();
                for (int k = 0; k < r.nextInt(4); ++k) hb[r.nextInt(hb.length)] = ACGT[r.nextInt(4)];
                if (r.nextInt(3) == 0) { // indel
                    int at = r.nextInt(hb.length), len = 1 + r.nextInt(8);
                    byte[] nb = new byte[hb.length + (r.nextBoolean() ? len : -Math.min(len, hb.length - at - 1))];
                    if (nb.length > hb.length) {
                        System.arraycopy(hb, 0, nb, 0, at);
                        for (int k = 0; k < len; ++k) nb[at + k] = ACGT[r.nextInt(4)];
                        System.arraycopy(hb, at, nb, at + len, hb.length - at);
                    } else {
                        int del = hb.length - nb.length;
                        System.arraycopy(hb, 0, nb, 0, at);
                        System.arraycopy(hb, at + del, nb, at, hb.length - at - del);
                    }
                    hb = nb;
                }
                if (r.nextInt(20) == 0) hb[r.nextInt(hb.length)] = 'N';
                haps.add(hb);
            }
            List<byte[][]> reads = new ArrayList<>();
            int nr = 1 + r.nextInt(6);
            for (int k = 0; k < nr; ++k) {
                byte[] src = haps.get(r.nextInt(nh));
                int len = 1 + r.nextInt(Math.min(160, src.length));
                if (r.nextInt(10) == 0) len = Math.min(160, src.length + r.nextInt(20)); // longer than the haplotype
                byte[] rb = new byte[len], q = new byte[len], ins = new byte[len], del = new byte[len], gcp = new byte[len];
                int off = src.length > len ? r.nextInt(src.length - len + 1) : 0;
                for (int i = 0; i < len; ++i) {
                    rb[i] = (off + i < src.length && r.nextInt(50) != 0) ? src[off + i] : ACGT[r.nextInt(4)];
                    if (r.nextInt(200) == 0) rb[i] = 'N';
                    q[i] = (byte) (r.nextInt(8) == 0 ? 6 : 18 + r.nextInt(24));
                    ins[i] = (byte) (r.nextInt(4) == 0 ? 6 + r.nextInt(40) : 45);
                    del[i] = (byte) (r.nextInt(4) == 0 ? 6 + r.nextInt(40) : 45);
                    gcp[i] = 10;
                }
                reads.add(new byte[][] {rb, q, ins, del, gcp});
            }
            int maxRead = 0, maxHap = 0;
            for (byte[][] rd : reads) maxRead = Math.max(maxRead, rd[0].length);
            for (byte[] h : haps) maxHap = Math.max(maxHap, h.length);
            LoglessPairHMM hmm = new LoglessPairHMM();
            hmm.initialize(maxRead, maxHap);
            for (byte[][] rd : reads) {
                for (int h = 0; h < nh; ++h) {
                    byte[] next = h + 1 < nh ? haps.get(h + 1) : null;
                    double lk = hmm.computeReadLikelihoodGivenHaplotypeLog10(haps.get(h), rd[0], rd[1], rd[2], rd[3], rd[4], true, next);
                    out.append(new String(haps.get(h))).append(' ').append(new String(rd[0])).append(' ').append(join(rd[1])).append(' ')
                       .append(join(rd[2])).append(' ').append(join(rd[3])).append(' ').append(join(rd[4])).append(' ')
                       .append(Long.toHexString(Double.doubleToRawLongBits(lk))).append('\n');
                    ++done;
                }
            }
        }
        System.out.print(out);
    }
}
