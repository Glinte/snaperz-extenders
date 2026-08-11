/* Round-18/21 audit: the canonical Xi-section of A_L, its return-weight word,
 * the weight palindrome, and the L=69 "UDDDD" local section criterion.
 *
 * The section eta_{L-2}(C_{L-2}) is characterised (verified separately, in
 * r18_eco.py) as: rightmost nonzero coordinate == 1, at least one zero after
 * it, and left neighbour >= 2.
 *
 * usage: r18_xi L [--dist]
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

static inline int lastpos(const uint8_t *x) {
    int i = L - 1;
    while (i >= 0 && !x[i]) i--;
    return i;
}

static inline int in_section(const uint8_t *x) {
    int j = lastpos(x);
    return j >= 1 && j + 1 < L && x[j] == 1 && x[j - 1] >= 2;
}

static int cmpl(const void *a, const void *b) {
    long x = *(const long *)a, y = *(const long *)b;
    return x < y ? -1 : x > y;
}

int main(int argc, char **argv) {
    L = atoi(argv[1]);
    int dist = argc > 2 && !strcmp(argv[2], "--dist");

    uint8_t *x = calloc(L, 1), *s0 = calloc(L, 1);
    /* anchor: A_L(R_L) = (L-1, 1, 0, ..., 0) */
    s0[0] = (uint8_t)(L - 1); s0[1] = 1;
    memcpy(x, s0, L);
    if (!in_section(x)) { printf("L=%d: anchor not in section\n", L); return 1; }

    long cap = 1 << 22, p = 0, nw = 0, last = 0;
    long *w = malloc(cap * sizeof(long));
    uint8_t *word = malloc(cap);
    int h = L - lastpos(x);
    for (;;) {
        sweep(x);
        p++;
        int h2 = L - lastpos(x);
        if (p > cap) { cap *= 2; w = realloc(w, cap * sizeof(long)); word = realloc(word, cap); }
        word[p - 1] = (uint8_t)(h2 - h + 1);   /* 0=D 1=F 2=U */
        h = h2;
        if (in_section(x)) { w[nw++] = p - last; last = p; }
        if (!memcmp(x, s0, L)) break;
    }

    int pal = 1;
    long sum = 0;
    for (long i = 0; i < nw; i++) {
        sum += w[i];
        if (w[i] != w[nw - 1 - i]) pal = 0;
    }
    printf("L=%d p=%ld Xi_period=%ld sum_w=%ld palindrome=%d\n",
           L, p, nw, sum, pal);

    if (dist) {
        long *c = malloc(nw * sizeof(long));
        memcpy(c, w, nw * sizeof(long));
        qsort(c, nw, sizeof(long), cmpl);
        printf("  weight distribution:");
        for (long i = 0; i < nw;) {
            long j = i;
            while (j < nw && c[j] == c[i]) j++;
            printf(" %ld^%ld", c[i], j - i);
            i = j;
        }
        printf("\n");
        /* UDDDD criterion: symbols at t-2..t+2 around a section time */
        uint8_t *sec = calloc(p, 1);
        long t = 0;
        for (long i = 0; i < nw; i++) { t += w[i]; sec[t % p] = 1; }
        long fp = 0, fn = 0, hits = 0;
        for (long u = 0; u < p; u++) {
            int ud = 1;
            for (int k = -2; k <= 2; k++) {
                long s = ((u + k) % p + p) % p;
                int want = (k == -2) ? 2 : 0;
                if (word[s] != want) { ud = 0; break; }
            }
            if (ud) hits++;
            if (ud && !sec[u]) fp++;
            if (!ud && sec[u]) fn++;
        }
        printf("  UDDDD windows=%ld  Xi section times=%ld  false_pos=%ld false_neg=%ld\n",
               hits, nw, fp, fn);
    }
    return 0;
}
