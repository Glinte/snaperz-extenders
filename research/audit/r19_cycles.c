/* Round-19 audit: complete cycle decomposition of A_n by direct enumeration,
 * to check the recursively reconstructed censuses for A_14 and A_15.
 *
 * usage: r19_cycles n
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n;
static uint64_t *states;
static long ns, cap;
static uint8_t cur[20];

static uint64_t pack(const uint8_t *x) {
    uint64_t v = 0;
    for (int i = 0; i < n; i++) v |= (uint64_t)x[i] << (4 * i);
    return v;
}

static void gen(int pos, int s) {
    if (pos == n) {
        if (s == n) {
            if (ns == cap) { cap *= 2; states = realloc(states, cap * 8); }
            states[ns++] = pack(cur);
        }
        return;
    }
    int lo = pos + 1 - s; if (lo < 0) lo = 0;
    for (int v = lo; v <= n - s; v++) { cur[pos] = (uint8_t)v; gen(pos + 1, s + v); }
    cur[pos] = 0;
}

static int cmp(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : x > y;
}

static long find(uint64_t v) {
    long lo = 0, hi = ns - 1;
    while (lo <= hi) {
        long mid = (lo + hi) / 2;
        if (states[mid] == v) return mid;
        if (states[mid] < v) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

int main(int argc, char **argv) {
    n = atoi(argv[1]);
    cap = 1 << 20; states = malloc(cap * 8);
    gen(0, 0);
    qsort(states, ns, 8, cmp);
    printf("n=%d |C_n|=%ld\n", n, ns);

    int32_t *succ = malloc(ns * sizeof(int32_t));
    for (long i = 0; i < ns; i++) {
        uint8_t x[20];
        uint64_t v = states[i];
        for (int j = 0; j < n; j++) x[j] = (v >> (4 * j)) & 15;
        for (int j = 0; j < n - 1; j++) {
            uint8_t a = x[j];
            if (!a) continue;
            if (a == 1) { x[j] = 1 + x[j + 1]; x[j + 1] = 0; }
            else { x[j] = a - 1; x[j + 1] += 1; }
        }
        long k = find(pack(x));
        if (k < 0) { printf("A left the state space!\n"); return 1; }
        succ[i] = (int32_t)k;
    }

    uint8_t *seen = calloc(ns, 1);
    long ncyc = 0, maxp = 0, minp = 1L << 60;
    long *lens = malloc(ns * sizeof(long)); long nl = 0;
    for (long i = 0; i < ns; i++) {
        if (seen[i]) continue;
        long len = 0, j = i;
        do { seen[j] = 1; j = succ[j]; len++; } while (j != i);
        ncyc++; lens[nl++] = len;
        if (len > maxp) maxp = len;
        if (len < minp) minp = len;
    }
    /* distinct period lengths */
    qsort(lens, nl, sizeof(long), (int (*)(const void *, const void *))cmp);
    long distinct = 0;
    for (long i = 0; i < nl; i++) if (!i || lens[i] != lens[i - 1]) distinct++;
    printf("cycles=%ld distinct_periods=%ld min=%ld max=%ld\n",
           ncyc, distinct, minp, maxp);
    return 0;
}
