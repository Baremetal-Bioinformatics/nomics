/* Write (x, erfc(x)) bit pairs from the system libm for nomics' erfc port.
 * usage: erfc_samples N > out.bin   (cc -O2 erfc_samples.c -lm) */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static uint64_t s = 0x9e3779b97f4a7c15ULL;
static uint64_t next(void) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s; }
int main(int argc, char **argv) {
	long n = atol(argv[1]), i;
	for (i = 0; i < n; ++i) {
		double x;
		uint64_t r = next(), b;
		switch (i % 4) {
		case 0: x = (double)(r >> 11) / 9007199254740992.0 * 6.0; break;        /* [0,6): bwa's range */
		case 1: x = (double)(r >> 11) / 9007199254740992.0 * 60.0 - 30.0; break; /* both signs, all branches */
		case 2: b = r & 0x800fffffffffffffULL; b |= (uint64_t)(0x3c0 + (next() % 0x50)) << 52; memcpy(&x, &b, 8); break; /* many exponents */
		default: x = fabs((double)((int64_t)(r % 4000) - 2000) / 37.3) * M_SQRT1_2; /* ns-like */
		}
		double y = erfc(x);
		fwrite(&x, 8, 1, stdout); fwrite(&y, 8, 1, stdout);
	}
	return 0;
}
