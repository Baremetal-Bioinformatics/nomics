#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../ksort.h"
#include "../kbtree.h"

typedef struct {
    int64_t key;
    uint32_t id;
    uint8_t payload[28];
} item_t;

_Static_assert(sizeof(item_t) == 40, "the BWA mem_chain_t oracle is 40 bytes");

#define item_lt(a, b) ((a).key < (b).key)
KSORT_INIT(item, item_t, item_lt)
#define item_gt(a, b) ((a).key > (b).key)
KSORT_INIT(itemdesc, item_t, item_gt)

static int item_cmp(item_t a, item_t b)
{
    return (a.key > b.key) - (a.key < b.key);
}
KBTREE_INIT(item, item_t, item_cmp)

static uint64_t next_u64(uint64_t *x)
{
    *x ^= *x << 13;
    *x ^= *x >> 7;
    *x ^= *x << 17;
    return *x;
}

/* The body is `ks_introsort_item` with this one counter added at its
 * ks_combsort call. It establishes that the adversarial case reaches klib's
 * depth-limit fallback; normal comparisons still use the header instantiation. */
static int counted_introsort(item_t a[], size_t n)
{
    int d, falls = 0;
    ks_isort_stack_t *top, *stack;
    item_t rp, swap_tmp;
    item_t *s, *t, *i, *j, *k;
    if (n < 1) return 0;
    else if (n == 2) { if (item_lt(a[1], a[0])) { swap_tmp = a[0]; a[0] = a[1]; a[1] = swap_tmp; } return 0; }
    for (d = 2; 1ul << d < n; ++d);
    stack = (ks_isort_stack_t*)malloc(sizeof(ks_isort_stack_t) * ((sizeof(size_t)*d)+2));
    top = stack; s = a; t = a + (n-1); d <<= 1;
    while (1) {
        if (s < t) {
            if (--d == 0) { ++falls; ks_combsort_item(t - s + 1, s); t = s; continue; }
            i = s; j = t; k = i + ((j-i)>>1) + 1;
            if (item_lt(*k, *i)) { if (item_lt(*k, *j)) k = j; }
            else k = item_lt(*j, *i)? i : j;
            rp = *k;
            if (k != t) { swap_tmp = *k; *k = *t; *t = swap_tmp; }
            for (;;) {
                do ++i; while (item_lt(*i, rp));
                do --j; while (i <= j && item_lt(rp, *j));
                if (j <= i) break;
                swap_tmp = *i; *i = *j; *j = swap_tmp;
            }
            swap_tmp = *i; *i = *t; *t = swap_tmp;
            if (i-s > t-i) {
                if (i-s > 16) { top->left = s; top->right = i-1; top->depth = d; ++top; }
                s = t-i > 16? i+1 : t;
            } else {
                if (t-i > 16) { top->left = i+1; top->right = t; top->depth = d; ++top; }
                t = i-s > 16? i-1 : s;
            }
        } else {
            if (top == stack) { free(stack); __ks_insertsort_item(a, a+n); return falls; }
            else { --top; s = (item_t*)top->left; t = (item_t*)top->right; d = top->depth; }
        }
    }
}

static void run_sort(int descending, uint64_t seed, size_t n, int kind)
{
    size_t i;
    item_t *a = (item_t*)calloc(n ? n : 1, sizeof(item_t));
    for (i = 0; i < n; ++i) {
        uint64_t r = next_u64(&seed);
        a[i].id = (uint32_t)i;
        if (kind == 0) a[i].key = (int64_t)(r % 7) - 3;
        else if (kind == 1) a[i].key = (int64_t)r;
        else a[i].key = (int64_t)((i & 1)? (n/2 + i/2) : i/2);
    }
    if (descending) ks_introsort(itemdesc, n, a);
    else ks_introsort(item, n, a);
    for (i = 0; i < n; ++i) printf("%u%c", a[i].id, i + 1 == n ? '\n' : ' ');
    if (!n) puts("");
    free(a);
}

static void run_fallback(size_t n)
{
    size_t i;
    item_t *a = (item_t*)calloc(n, sizeof(item_t));
    for (i = 0; i < n; ++i) {
        /* Median-of-three killer: odds followed by evens. */
        a[i].key = (int64_t)((i & 1)? n/2 + i/2 : i/2);
        a[i].id = (uint32_t)i;
    }
    printf("%d\n", counted_introsort(a, n));
    free(a);
}

static void visit(item_t *x) { printf("%u ", x->id); }
static void run_btree(int size, uint64_t seed, int ops)
{
    int i;
    kbtree_t(item) *b = kb_init(item, size);
    for (i = 0; i < ops; ++i) {
        item_t query, value, *lower, *upper;
        memset(&query, 0, sizeof(query));
        query.key = (int64_t)(next_u64(&seed) % 23) - 11;
        kb_intervalp(item, b, &query, &lower, &upper);
        printf("Q %d %d\n", lower ? (int)lower->id : -1, upper ? (int)upper->id : -1);
        if (lower && ((lower->id ^ (uint32_t)next_u64(&seed)) & 1)) {
            lower->payload[0] ^= 0x5a;
        } else {
            memset(&value, 0, sizeof(value));
            value.key = query.key; value.id = (uint32_t)i;
            kb_putp(item, b, &value);
        }
    }
    printf("T ");
    __kb_traverse(item_t, b, visit);
    puts("");
    kb_destroy(item, b);
}

int main(int argc, char **argv)
{
    if (argc < 2) return 2;
    if (!strcmp(argv[1], "sort") && argc == 6) {
        run_sort(atoi(argv[2]), strtoull(argv[3], 0, 10), strtoull(argv[4], 0, 10), atoi(argv[5]));
    } else if (!strcmp(argv[1], "fallback") && argc == 3) {
        run_fallback(strtoull(argv[2], 0, 10));
    } else if (!strcmp(argv[1], "btree") && argc == 5) {
        run_btree(atoi(argv[2]), strtoull(argv[3], 0, 10), atoi(argv[4]));
    } else return 2;
    return 0;
}
