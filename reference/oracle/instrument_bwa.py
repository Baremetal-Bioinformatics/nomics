#!/usr/bin/env python3
"""Add stage dumps to a pristine bwa v0.7.19 source tree (in place).

usage: instrument_bwa.py BWA_SRC_DIR

With BWA_DUMP=<dir> set, each worker thread appends lines to
<dir>/dump.<thread>.txt. Every line starts with a tag and the read id
`pair_index << 1 | mate` (pair_index counts from 0 over the whole run):

  M id  qb,qe,x0,x1,x2 ...           SMEM/seed intervals after mem_collect_intv
  c id  ci pos rid n | rbeg,qbeg,len,score ...   chains from mem_chain (kbtree order)
  C id  ci pos rid kept first w frac_rep n | seeds  chains after mem_chain_flt
  R id  ri rid rb re qb qe score truesc sub alt_sc csub sub_n w seedcov
        secondary seedlen0 n_comp frac_rep   regions after mem_sort_dedup_patch

Floats are printed with %a (exact). Without BWA_DUMP the binary behaves as
the pristine one. Test oracle only (docs/GATES.md).
"""
import pathlib
import sys

src = pathlib.Path(sys.argv[1]) / "bwamem.c"
s = src.read_text()


def patch(old, new, count=1):
    global s
    assert s.count(old) == count, (old, s.count(old))
    s = s.replace(old, new)


HELPERS = r'''
/* ---- nomics oracle dump ---- */
static __thread FILE *nd_fp;
__thread long nd_id = -1;
static FILE *nd_file(void)
{
	if (!nd_fp) {
		const char *d = getenv("BWA_DUMP");
		char p[4096];
		if (!d) return 0;
		snprintf(p, sizeof p, "%s/dump.%lx.txt", d, (unsigned long)pthread_self());
		nd_fp = fopen(p, "a");
	}
	return nd_fp;
}
'''

patch('#include "kbtree.h"\n', '#include "kbtree.h"\n' + HELPERS)

# SMEMs
patch(
    "\tmem_collect_intv(opt, bwt, len, seq, aux);\n",
    "\tmem_collect_intv(opt, bwt, len, seq, aux);\n"
    "\tif (nd_file()) { FILE *f = nd_file(); size_t j; fprintf(f, \"M %ld\", nd_id);\n"
    "\t\tfor (j = 0; j < aux->mem.n; ++j) { bwtintv_t *q = &aux->mem.a[j];\n"
    "\t\t\tfprintf(f, \" %d,%d,%ld,%ld,%ld\", (int)(q->info>>32), (int)(uint32_t)q->info, (long)q->x[0], (long)q->x[1], (long)q->x[2]); }\n"
    "\t\tfputc('\\n', f); }\n",
)

CHAIN_DUMP = r'''
static void nd_chains(const char *tag, const mem_chain_v *chn, int full)
{
	FILE *f = nd_file();
	size_t i;
	int j;
	if (!f) return;
	for (i = 0; i < chn->n; ++i) {
		const mem_chain_t *c = &chn->a[i];
		if (full) fprintf(f, "%s %ld %zu %ld %d %d %d %d %a %d |", tag, nd_id, i, (long)c->pos, c->rid, (int)c->kept, c->first, (int)c->w, c->frac_rep, c->n);
		else fprintf(f, "%s %ld %zu %ld %d %d |", tag, nd_id, i, (long)c->pos, c->rid, c->n);
		for (j = 0; j < c->n; ++j) fprintf(f, " %ld,%d,%d,%d", (long)c->seeds[j].rbeg, c->seeds[j].qbeg, c->seeds[j].len, c->seeds[j].score);
		fputc('\n', f);
	}
}

static void nd_regs(const mem_alnreg_v *r)
{
	FILE *f = nd_file();
	size_t i;
	if (!f) return;
	for (i = 0; i < r->n; ++i) {
		const mem_alnreg_t *p = &r->a[i];
		fprintf(f, "R %ld %zu %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %a\n", nd_id, i, p->rid, (long)p->rb, (long)p->re, p->qb, p->qe,
			p->score, p->truesc, p->sub, p->alt_sc, p->csub, p->sub_n, p->w, p->seedcov, p->secondary, p->seedlen0, (int)p->n_comp, p->frac_rep);
	}
	fflush(f);
}

mem_alnreg_v mem_align1_core('''

patch("\nmem_alnreg_v mem_align1_core(", CHAIN_DUMP)
patch(
    "\tchn = mem_chain(opt, bwt, bns, l_seq, (uint8_t*)seq, buf);\n\tchn.n = mem_chain_flt(opt, chn.n, chn.a);\n",
    "\tchn = mem_chain(opt, bwt, bns, l_seq, (uint8_t*)seq, buf);\n\tnd_chains(\"c\", &chn, 0);\n"
    "\tchn.n = mem_chain_flt(opt, chn.n, chn.a);\n",
)
patch(
    "\tif (bwa_verbose >= 4) mem_print_chain(bns, &chn);\n\n\tkv_init(regs);",
    "\tif (bwa_verbose >= 4) mem_print_chain(bns, &chn);\n\tnd_chains(\"C\", &chn, 1);\n\n\tkv_init(regs);",
)
patch(
    "\t\t\tp->is_alt = 1;\n\t}\n\treturn regs;\n}",
    "\t\t\tp->is_alt = 1;\n\t}\n\tnd_regs(&regs);\n\treturn regs;\n}",
)

# Read ids in worker1 (paired-end path only).
patch(
    "\t\tw->regs[i<<1|0] = mem_align1_core(",
    "\t\tnd_id = (long)(((w->n_processed>>1) + i) << 1 | 0);\n\t\tw->regs[i<<1|0] = mem_align1_core(",
)
patch(
    "\t\tw->regs[i<<1|1] = mem_align1_core(",
    "\t\tnd_id = (long)(((w->n_processed>>1) + i) << 1 | 1);\n\t\tw->regs[i<<1|1] = mem_align1_core(",
)

src.write_text(s)
print("instrumented", src)
