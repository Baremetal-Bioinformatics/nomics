/* Batch protocol for ksw_align2_vs_c.rs. ksw.c itself is unmodified. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "ksw.h"

static int32_t get_i32(FILE *f)
{
	int32_t x;
	if (fread(&x, sizeof(x), 1, f) != 1) exit(2);
	return x;
}

static void put_i32(FILE *f, int32_t x)
{
	if (fwrite(&x, sizeof(x), 1, f) != 1) exit(3);
}

int main(int argc, char **argv)
{
	FILE *in, *out;
	uint32_t n, case_no;
	if (argc != 3) return 64;
	in = fopen(argv[1], "rb");
	out = fopen(argv[2], "wb");
	if (!in || !out || fread(&n, sizeof(n), 1, in) != 1) return 65;
	for (case_no = 0; case_no < n; ++case_no) {
		int32_t qlen = get_i32(in), tlen = get_i32(in), m = get_i32(in);
		int32_t o_del = get_i32(in), e_del = get_i32(in);
		int32_t o_ins = get_i32(in), e_ins = get_i32(in), xtra = get_i32(in);
		int8_t *mat = malloc((size_t)m * m);
		uint8_t *query = malloc((size_t)qlen);
		uint8_t *target = malloc((size_t)tlen);
		kswr_t r;
		if (!mat || !query || (tlen && !target)) return 66;
		if (fread(mat, 1, (size_t)m * m, in) != (size_t)m * m ||
		    fread(query, 1, (size_t)qlen, in) != (size_t)qlen ||
		    fread(target, 1, (size_t)tlen, in) != (size_t)tlen) return 67;
		r = ksw_align2(qlen, query, tlen, target, m, mat,
		               o_del, e_del, o_ins, e_ins, xtra, 0);
		put_i32(out, r.score); put_i32(out, r.te); put_i32(out, r.qe);
		put_i32(out, r.score2); put_i32(out, r.te2);
		put_i32(out, r.tb); put_i32(out, r.qb);
		free(mat); free(query); free(target);
	}
	if (fclose(out) || fclose(in)) return 68;
	return 0;
}
