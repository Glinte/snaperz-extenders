/* Round-20 audit: the final-descent (U/F/D) itinerary of the canonical orbit,
 * its exact LZ78 phrase count, and the claimed global height-reflection law.
 *
 * usage: r20_lz L [--word out.bin]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int L;

static inline void sweep(uint8_t *x) {
    for (int i = 0; i < L - 1; i++) {
        uint8_t a = x[i];
        if (!a) continue;
        if (a == 1) { x[i] = 1 + x[i + 1]; x[i + 1] = 0; }
        else { x[i] = a - 1; x[i + 1] += 1; }
    }
}

static inline int height(const uint8_t *x) {
    int i = L - 1;
    while (i >= 0 && x[i] == 0) i--;
    return L - i;                 /* length of the final descent */
}

/* LZ78: count phrases of a ternary word. */
static long lz78_phrases(const uint8_t *w, long n) {
    long cap = 1 << 20, cnt = 1;   /* node 0 = root */
    int32_t *child = malloc((size_t)cap * 3 * sizeof(int32_t));
    memset(child, 0, (size_t)cap * 3 * sizeof(int32_t));
    long phrases = 0, cur = 0;
    for (long i = 0; i < n; i++) {
        int c = w[i];
        int32_t nx = child[cur * 3 + c];
        if (nx) { cur = nx; continue; }
        if (cnt == cap) {
            long ncap = cap * 2;
            child = realloc(child, (size_t)ncap * 3 * sizeof(int32_t));
            memset(child + cap * 3, 0, (size_t)cap * 3 * sizeof(int32_t));
            cap = ncap;
        }
        child[cur * 3 + c] = (int32_t)cnt++;
        phrases++;
        cur = 0;
    }
    if (cur) phrases++;            /* trailing incomplete phrase */
    free(child);
    return phrases;
}

int main(int argc, char **argv) {
    L = atoi(argv[1]);
    const char *wordfile = NULL;
    for (int i = 2; i < argc - 1; i++)
        if (!strcmp(argv[i], "--word")) wordfile = argv[i + 1];

    uint8_t *x = malloc(L), *x0 = malloc(L);
    for (int i = 0; i < L; i++) x0[i] = 1;
    memcpy(x, x0, L);

    long cap = 1 << 22, p = 0;
    uint8_t *w = malloc(cap);      /* 0=D 1=F 2=U */
    int h = height(x);
    for (;;) {
        sweep(x);
        p++;
        int h2 = height(x);
        if (p > cap) { cap *= 2; w = realloc(w, cap); }
        w[p - 1] = (uint8_t)(h2 - h + 1);
        h = h2;
        if (!memcmp(x, x0, L)) break;
    }

    /* height reflection: delta_t = -delta_{p - 2L + 1 - t mod p} */
    long bad = 0, off = ((p - 2 * (long)L + 1) % p + p) % p;
    for (long t = 0; t < p; t++) {
        long s = ((off - t) % p + p) % p;
        if (w[t] + w[s] != 2) bad++;
    }

    long g = lz78_phrases(w, p);
    printf("L=%d p=%ld lz78=%ld frac=%.4f%% reflection_failures=%ld\n",
           L, p, g, 100.0 * g / p, bad);
    if (wordfile) { FILE *f = fopen(wordfile, "wb"); fwrite(w, 1, p, f); fclose(f); }
    return 0;
}
