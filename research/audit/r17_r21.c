/* Round-17 / Round-21 audit.
 *
 *  --defect L : the odd "defect gas" statistic delta(x) over the canonical
 *               orbit (claim: max delta = (L-15)/2 for odd 17 <= L <= 45).
 *  --tiles    : at L=69, the sweeps whose cell 55 is zero and whose 13-cell
 *               tail sits on the canonical A_13 cycle (claim: 514 maximal runs
 *               of 1419 = 11*129 sweeps, 729,366 sweeps, 64.8319%).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int L;

static inline void sweep_n(uint8_t *x, int n) {
    for (int i = 0; i < n - 1; i++) {
        uint8_t a = x[i];
        if (!a) continue;
        if (a == 1) { x[i] = 1 + x[i + 1]; x[i + 1] = 0; }
        else { x[i] = a - 1; x[i + 1] += 1; }
    }
}

static inline int good(int u, int v) {
    if (u == 0) return v % 2 == 0;
    return (u % 2 == 1) && (v % 2 == 1);
}

static void defect(int Larg) {
    L = Larg;
    uint8_t *x = calloc(L, 1), *x0 = calloc(L, 1);
    for (int i = 0; i < L; i++) x0[i] = 1;
    memcpy(x, x0, L);
    long p = 0, best = 0, zero = 0;
    double sum = 0;
    do {
        int mn = 1 << 30;
        for (int k = 0; k < L; k += 2) {
            int bad = 0;
            for (int i = 0; i + 1 < k; i += 2) if (!good(x[i], x[i + 1])) bad++;
            for (int i = k + 1; i + 1 < L; i += 2) if (!good(x[i], x[i + 1])) bad++;
            if (bad < mn) mn = bad;
        }
        if (mn > best) best = mn;
        if (mn == 0) zero++;
        sum += mn;
        sweep_n(x, L);
        p++;
    } while (memcmp(x, x0, L));
    printf("L=%-3d p=%-9ld max_delta=%-3ld  (L-15)/2=%-3d  mean=%.4f  zero=%.4f%%\n",
           L, p, best, (L - 15) / 2, sum / p, 100.0 * zero / p);
    free(x); free(x0);
}

static void tiles(void) {
    L = 69;
    const int S = 13, CUT = 56;
    uint8_t cyc[129][13];
    uint8_t t[13];
    for (int i = 0; i < S; i++) t[i] = 1;
    for (int k = 0; k < 129; k++) { memcpy(cyc[k], t, S); sweep_n(t, S); }

    uint8_t *x = calloc(L, 1), *x0 = calloc(L, 1);
    for (int i = 0; i < L; i++) x0[i] = 1;
    memcpy(x, x0, L);
    long p = 0, hit = 0, runs = 0, curlen = 0, minrun = 1 << 30, maxrun = 0;
    long len1419 = 0;
    int prev = 0;
    do {
        int ok = (x[CUT - 1] == 0);
        if (ok) {
            int f = 0;
            for (int k = 0; k < 129 && !f; k++)
                if (!memcmp(cyc[k], x + CUT, S)) f = 1;
            ok = f;
        }
        if (ok) { hit++; curlen++; }
        else if (prev) {
            runs++;
            if (curlen < minrun) minrun = curlen;
            if (curlen > maxrun) maxrun = curlen;
            if (curlen == 1419) len1419++;
            curlen = 0;
        }
        prev = ok;
        sweep_n(x, L);
        p++;
    } while (memcmp(x, x0, L));
    if (prev) { runs++; if (curlen == 1419) len1419++; if (curlen > maxrun) maxrun = curlen; }
    printf("L=69 p=%ld  synchronised sweeps=%ld (%.4f%%)  runs=%ld  "
           "runs_of_1419=%ld  min=%ld max=%ld\n",
           p, hit, 100.0 * hit / p, runs, len1419, minrun, maxrun);
}

int main(int argc, char **argv) {
    if (argc > 1 && !strcmp(argv[1], "--tiles")) { tiles(); return 0; }
    if (argc > 2 && !strcmp(argv[1], "--defect")) { defect(atoi(argv[2])); return 0; }
    for (int l = 17; l <= 45; l += 2) defect(l);
    return 0;
}
