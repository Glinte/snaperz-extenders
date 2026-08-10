/* Round-15 audit, part B: the general phase-reflection constant.
 *
 * Round 14 logged, at L=69 with the clean 56|13 cut:
 *     phi_j + phi_-j = z(P_j) + 43  (mod 129).
 * The round-15 write-up proposes the general law
 *     phi_j + phi_-j = z(P_j) + n + B_s - 1  (mod p_s)
 * for any clean n|s cut.  Since p_s = B_s + s - 1 this is just
 *     const = n - s  (mod p_s),
 * with B_s cancelling.  This program measures the actual constant over every
 * clean cut of every orbit it can run, and compares.
 *
 * Build: cc -O2 -o r15b r15b.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXL 80

static inline void Astep(unsigned char *x, int L) {
    for (int i = 0; i < L - 1; i++) {
        int a = x[i];
        if (!a) continue;
        if (a == 1) { x[i] = (unsigned char)(1 + x[i + 1]); x[i + 1] = 0; }
        else { x[i] = (unsigned char)(a - 1); x[i + 1]++; }
    }
}

/* ---- canonical A_s cycle, with a phase lookup ---- */
static unsigned char *cycbuf;
static int cyclen, cyc_s;
static int *htab, hmask;

static unsigned long hashv(const unsigned char *t, int s) {
    unsigned long h = 1469598103934665603UL;
    for (int i = 0; i < s; i++) { h ^= t[i]; h *= 1099511628211UL; }
    return h;
}

static void build_cycle(int s) {
    unsigned char y[MAXL], y0[MAXL];
    int cap = 1 << 12;
    cycbuf = malloc((size_t)cap * s);
    for (int i = 0; i < s; i++) y[i] = y0[i] = 1;
    cyclen = 0;
    for (;;) {
        if (cyclen == cap) { cap *= 2; cycbuf = realloc(cycbuf, (size_t)cap * s); }
        memcpy(cycbuf + (size_t)cyclen * s, y, s);
        cyclen++;
        Astep(y, s);
        if (!memcmp(y, y0, s)) break;
    }
    cyc_s = s;
    int hb = 1; while ((1 << hb) < cyclen * 4) hb++;
    hmask = (1 << hb) - 1;
    htab = malloc(sizeof(int) * (hmask + 1));
    for (int i = 0; i <= hmask; i++) htab[i] = -1;
    for (int i = 0; i < cyclen; i++) {
        unsigned long h = hashv(cycbuf + (size_t)i * s, s) & hmask;
        while (htab[h] >= 0) h = (h + 1) & hmask;
        htab[h] = i;
    }
}

static int phase_of(const unsigned char *t) {
    int s = cyc_s;
    unsigned long h = hashv(t, s) & hmask;
    while (htab[h] >= 0) {
        if (!memcmp(cycbuf + (size_t)htab[h] * s, t, s)) return htab[h];
        h = (h + 1) & hmask;
    }
    return -1;
}

static void free_cycle(void) { free(cycbuf); free(htab); }

/* ---- checkpoint table of the canonical orbit ---- */
static unsigned char *cps;   /* (N+1) x L, cps[0] = cps[N] = E_L */
static long Ncp;

static long build_orbit(int L) {
    unsigned char x[MAXL], x0[MAXL];
    for (int i = 0; i < L; i++) x[i] = x0[i] = 1;
    long cap = 1 << 12, k = 0, p = 0;
    cps = malloc((size_t)cap * L);
    memcpy(cps, x, L);            /* k = 0 : E_L */
    for (;;) {
        Astep(x, L); p++;
        if (x[0] == 1) {
            k++;
            if (k == cap) { cap *= 2; cps = realloc(cps, (size_t)cap * L); }
            memcpy(cps + (size_t)k * L, x, L);
        }
        if (!memcmp(x, x0, L)) break;
    }
    Ncp = k;                      /* cps[N] = E_L again */
    return p;
}

static int zdepth(const unsigned char *x, int n) {
    for (int i = 0; i < n; i++) if (!x[i]) return i;
    return n;
}

/* test the law for one clean cut n|s of one L; block = 1 to restrict to
 * 112-checkpoint sections (round 14's setup), 0 to use every checkpoint. */
static void test_cut(int L, int n, int block) {
    int s = L - n;
    build_cycle(s);
    int ps = cyclen;
    long N = Ncp;

    char *clean = calloc(N, 1);
    int *phi = malloc(sizeof(int) * N), *zp = malloc(sizeof(int) * N);
    long nclean = 0, offcyc = 0, ncand = 0;

    for (long k = 0; k < N; k++) {
        if (block && (k % 112)) continue;
        const unsigned char *x = cps + (size_t)k * L;
        int sum = 0;
        for (int i = 0; i < n; i++) sum += x[i];
        if (sum != n) continue;
        ncand++;
        int ph = phase_of(x + n);
        if (ph < 0) { offcyc++; continue; }
        clean[k] = 1; nclean++;
        phi[k] = ph; zp[k] = zdepth(x, n);
    }

    long pairs = 0, zfail = 0, cfail = 0;
    int c0 = -1;
    int cpred = ((n - s) % ps + ps) % ps;
    for (long k = 1; k < N; k++) {
        long m = (N - k) % N;
        if (!clean[k] || !clean[m]) continue;
        pairs++;
        if (zp[m] != zp[k]) zfail++;
        int c = ((phi[k] + phi[m] - zp[k]) % ps + ps) % ps;
        if (c0 < 0) c0 = c; else if (c != c0) cfail++;
    }
    printf("  L=%2d  n=%2d s=%2d  p_s=%6d  clean cuts %6ld (off-cycle %5ld of %6ld)  "
           "pairs %6ld  const %s",
           L, n, s, ps, nclean, offcyc, ncand, pairs,
           pairs == 0 ? "n/a" : (cfail ? "VARIES" : "fixed"));
    if (pairs && !cfail)
        printf(" = %3d   predicted n-s = %3d  %s", c0, cpred,
               c0 == cpred ? "MATCH" : "*** MISMATCH ***");
    if (pairs) printf("   z(P_-j)=z(P_j) fails %ld", zfail);
    printf("\n");
    free(clean); free(phi); free(zp); free_cycle();
}

int main(int argc, char **argv) {
    int Ls[] = {17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 18, 20, 22, 24, 26,
                28, 30, 32, 34, 36, 38, 40, 42, 44, 69};
    int nL = (int)(sizeof(Ls) / sizeof(Ls[0]));
    int only69 = (argc > 1 && !strcmp(argv[1], "69"));

    if (argc == 3) {                       /* single cut: r15b L n */
        int L = atoi(argv[1]), n = atoi(argv[2]);
        build_orbit(L);
        test_cut(L, n, 0);
        free(cps);
        return 0;
    }

    for (int i = 0; i < nL; i++) {
        int L = Ls[i];
        if (only69 && L != 69) continue;
        long p = build_orbit(L);
        printf("L=%2d  p=%ld  B=%ld  checkpoints N=%ld  sections N/112=%ld\n",
               L, p, p - L + 1, Ncp, Ncp / 112);
        /* every cut that leaves a tail small enough to cycle cheaply */
        for (int s = 2; s <= 24 && s <= L - 2; s++) {
            int n = L - s;
            test_cut(L, n, 0);
        }
        if (L == 69) {
            printf("  -- restricted to the 112-block sections (round-14 setup) --\n");
            for (int s = 9; s <= 15; s += 2) test_cut(L, L - s, 1);
        }
        free(cps);
        fflush(stdout);
    }
    return 0;
}
