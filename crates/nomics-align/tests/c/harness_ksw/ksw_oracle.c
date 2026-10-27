/* Batch protocol for ksw_vs_c.rs.  ksw.c itself remains unmodified. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "ksw.h"

static int32_t get_i32(FILE *f) { int32_t x; if (fread(&x, sizeof(x), 1, f) != 1) exit(2); return x; }
static void put_i32(FILE *f, int32_t x) { if (fwrite(&x, sizeof(x), 1, f) != 1) exit(3); }
static void put_u32(FILE *f, uint32_t x) { if (fwrite(&x, sizeof(x), 1, f) != 1) exit(3); }

int main(int argc, char **argv)
{
	FILE *in, *out;
	uint32_t n, case_no;
	if (argc != 3) return 64;
	in = fopen(argv[1], "rb"); out = fopen(argv[2], "wb");
	if (!in || !out || fread(&n, sizeof(n), 1, in) != 1) return 65;
	for (case_no = 0; case_no < n; ++case_no) {
		uint8_t kind;
		int32_t qlen, tlen, m, o_del, e_del, o_ins, e_ins, w, end_bonus, zdrop, h0;
		uint8_t *query, *target; int8_t *mat;
		if (fread(&kind, 1, 1, in) != 1) return 66;
		qlen=get_i32(in); tlen=get_i32(in); m=get_i32(in); o_del=get_i32(in); e_del=get_i32(in);
		o_ins=get_i32(in); e_ins=get_i32(in); w=get_i32(in); end_bonus=get_i32(in); zdrop=get_i32(in); h0=get_i32(in);
		mat = malloc((size_t)m*m); query = malloc((size_t)qlen); target = malloc((size_t)tlen);
		if (fread(mat, 1, (size_t)m*m, in) != (size_t)m*m || fread(query, 1, qlen, in) != (size_t)qlen || fread(target, 1, tlen, in) != (size_t)tlen) return 67;
		if (kind == 0) {
			int qle, tle, gtle, gscore, max_off;
			int score = ksw_extend2(qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, w, end_bonus, zdrop, h0, &qle, &tle, &gtle, &gscore, &max_off);
			put_i32(out, score); put_i32(out, qle); put_i32(out, tle); put_i32(out, gtle); put_i32(out, gscore); put_i32(out, max_off);
		} else {
			int n_cigar, score, i; uint32_t *cigar;
			score = ksw_global2(qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, w, &n_cigar, &cigar);
			put_i32(out, score); put_u32(out, (uint32_t)n_cigar);
			for (i = 0; i < n_cigar; ++i) put_u32(out, cigar[i]);
			free(cigar);
		}
		free(mat); free(query); free(target);
	}
	if (fclose(out) || fclose(in)) return 68;
	return 0;
}
