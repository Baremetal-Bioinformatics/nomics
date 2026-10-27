/* Fixtures for nomics' glibc-exact log/exp: "<fn> <hex x> <hex result>" per line.
 * Built and run on the Spark (Ubuntu glibc 2.39, aarch64):
 *   cc -O2 glibc_logexp_samples.c -lm -o s && ./s N | gzip -9 > glibc_logexp.tsv.gz */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static uint64_t s = 0x243f6a8885a308d3ULL;
static uint64_t next(void) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s; }
static double unit(void) { return (double)(next() >> 11) / 9007199254740992.0; }
static uint64_t bits(double d) { uint64_t b; memcpy(&b, &d, 8); return b; }
static void out(const char *fn, double x, double y) { printf("%s %llx %llx\n", fn, (unsigned long long)bits(x), (unsigned long long)bits(y)); }
int main(int argc, char **argv) {
	long n = atol(argv[1]), i;
	for (i = 0; i < n; ++i) {
		double x;
		switch (i % 6) {
		case 0: x = unit() * 1000; break;                                  /* bwa: log(read length), log(n+1) */
		case 1: x = 1 + (unit() - 0.5) * 0.2; break;                       /* near 1: log's special path */
		case 2: { uint64_t b = next() & 0x7fefffffffffffffULL; memcpy(&x, &b, 8); break; } /* all positive finite */
		case 3: x = unit() * 2; break;                                     /* 2*erfc(..) in (0,2] */
		case 4: x = ldexp(unit() + 0.5, (int)(next() % 2100) - 1074); break; /* subnormal..huge */
		default: x = (double)(next() % 100000) + 1; break;                 /* integers */
		}
		out("log", x, log(x));
		double e;
		switch (i % 4) {
		case 0: e = (unit() - 0.5) * 1500; break;                          /* full range incl. over/underflow */
		case 1: e = -unit() * 50; break;                                   /* erfc's exp(-z*z-0.5625) */
		case 2: e = (unit() - 0.5) * 1e-6; break;                          /* tiny */
		default: { uint64_t b = next(); memcpy(&e, &b, 8); if (isnan(e)) e = 0.5; break; }
		}
		out("exp", e, exp(e));
	}
	double edge[] = {0.0, -0.0, 1.0, -1.0, INFINITY, -INFINITY, 4.9e-324, 2.2250738585072014e-308, 1.7976931348623157e308,
	                 709.782712893384, 709.7827128933841, -745.1332191019411, -745.1332191019412, 0.5, 2.0, 0x1.0p-1022, 0x1.fffffffffffffp-1};
	for (unsigned k = 0; k < sizeof edge / sizeof *edge; ++k) { out("log", edge[k], log(edge[k])); out("exp", edge[k], exp(edge[k])); }
	return 0;
}
