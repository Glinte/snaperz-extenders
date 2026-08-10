/* Round-15 audit, part C: where does the reflection constant actually come from?
 *
 * Round 14 proved, by telescoping, the checkpoint-time reflection identity
 *     t_k + t_{N-k} = p_L - L + z(x_k).
 * If on clean n|s cuts the tail phase is locked to the clock,
 *     phi_k = t_k + c   (mod p_s)   for a single offset c,
 * then immediately
 *     phi_k + phi_{N-k} = z(P_k) + (p_L - L + 2c)   (mod p_s),
 * i.e. the reflection constant is FORCED to be  p_L - L + 2c  mod p_s.
 * That would prove the law rather than fit it.  This measures c over every
 * clean checkpoint (round 14 only checked the 112-sections, j <= 780).
 *
 * Build: cc -O2 -o r15c r15c.c
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

static unsigned char *cycbuf; static int cyclen, cyc_s; static int *htab, hmask;
static unsigned long hashv(const unsigned char *t, int s) {
    unsigned long h = 1469598103934665603UL;
    for (int i = 0; i < s; i++) { h ^= t[i]; h *= 1099511628211UL; }
    return h;
}
static void build_cycle(int s) {
    unsigned char y[MAXL], y0[MAXL];
    int cap = 1 << 12; cycbuf = malloc((size_t)cap * s);
    for (int i = 0; i < s; i++) y[i] = y0[i] = 1;
    cyclen = 0;
    for (;;) {
        if (cyclen == cap) { cap *= 2; cycbuf = realloc(cycbuf, (size_t)cap * s); }
        memcpy(cycbuf + (size_t)cyclen * s, y, s); cyclen++;
        Astep(y, s); if (!memcmp(y, y0, s)) break;
    }
    cyc_s = s;
    int hb = 1; while ((1 << hb) < cyclen * 4) hb++;
    hmask = (1 << hb) - 1; htab = malloc(sizeof(int) * (hmask + 1));
    for (int i = 0; i <= hmask; i++) htab[i] = -1;
    for (int i = 0; i < cyclen; i++) {
        unsigned long h = hashv(cycbuf + (size_t)i * s, s) & hmask;
        while (htab[h] >= 0) h = (h + 1) & hmask;
        htab[h] = i;
    }
}
static int phase_of(const unsigned char *t) {
    unsigned long h = hashv(t, cyc_s) & hmask;
    while (htab[h] >= 0) {
        if (!memcmp(cycbuf + (size_t)htab[h] * cyc_s, t, cyc_s)) return htab[h];
        h = (h + 1) & hmask;
    }
    return -1;
}
static void free_cycle(void) { free(cycbuf); free(htab); }

static unsigned char *cps; static long *tks; static long Ncp, pL;

static void build_orbit(int L) {
    unsigned char x[MAXL], x0[MAXL];
    for (int i = 0; i < L; i++) x[i] = x0[i] = 1;
    long cap = 1 << 12, k = 0, p = 0;
    cps = malloc((size_t)cap * L); tks = malloc(sizeof(long) * cap);
    memcpy(cps, x, L); tks[0] = 0;
    for (;;) {
        Astep(x, L); p++;
        if (x[0] == 1) {
            k++;
            if (k == cap) { cap *= 2; cps = realloc(cps, (size_t)cap * L);
                            tks = realloc(tks, sizeof(long) * cap); }
            memcpy(cps + (size_t)k * L, x, L); tks[k] = p;
        }
        if (!memcmp(x, x0, L)) break;
    }
    Ncp = k; pL = p;
}

static int zdepth(const unsigned char *x, int n) {
    for (int i = 0; i < n; i++) if (!x[i]) return i;
    return n;
}

static void analyse(int L, int n) {
    int s = L - n;
    build_cycle(s);
    int ps = cyclen;
    long N = Ncp;
    char *clean = calloc(N, 1);
    int *phi = malloc(sizeof(int) * N), *zp = malloc(sizeof(int) * N);
    long nclean = 0, offcyc = 0;
    int *chist = calloc(ps, sizeof(int));

    for (long k = 0; k < N; k++) {
        const unsigned char *x = cps + (size_t)k * L;
        int sum = 0; for (int i = 0; i < n; i++) sum += x[i];
        if (sum != n) continue;
        int ph = phase_of(x + n);
        if (ph < 0) { offcyc++; continue; }
        clean[k] = 1; nclean++; phi[k] = ph; zp[k] = zdepth(x, n);
        chist[(int)((((phi[k] - tks[k]) % ps) + ps) % ps)]++;
    }
    if (nclean < 8) { free(clean); free(phi); free(zp); free(chist); free_cycle(); return; }

    int ndist = 0, cmode = -1, cmax = 0;
    for (int i = 0; i < ps; i++) if (chist[i]) { ndist++; if (chist[i] > cmax) { cmax = chist[i]; cmode = i; } }

    long pairs = 0, cfail = 0; int c0 = -1;
    for (long k = 1; k < N; k++) {
        long m = (N - k) % N;
        if (!clean[k] || !clean[m]) continue;
        pairs++;
        int c = ((phi[k] + phi[m] - zp[k]) % ps + ps) % ps;
        if (c0 < 0) c0 = c; else if (c != c0) cfail++;
    }
    int forced = (int)(((pL - L + 2L * cmode) % ps + ps) % ps);
    int naive  = (int)(((long)(n - s) % ps + ps) % ps);

    printf("  L=%2d n=%2d s=%2d p_s=%5d | clean %7ld off %6ld | lock offsets %4d distinct "
           "(mode %4d, %5.1f%%) | pairs %7ld const %s",
           L, n, s, ps, nclean, offcyc, ndist, cmode, 100.0 * cmax / (double)nclean,
           pairs, pairs == 0 ? "n/a" : cfail ? "VARIES" : "fixed");
    if (pairs && !cfail)
        printf(" %4d | p-L+2c = %4d %s | n-s = %4d %s", c0,
               forced, forced == c0 ? "OK " : "no ", naive, naive == c0 ? "OK " : "no ");
    printf("\n");
    free(clean); free(phi); free(zp); free(chist); free_cycle();
}

int main(void) {
    printf("lock offset c_k = (phi_k - t_k) mod p_s over EVERY clean checkpoint\n\n");
    int Ls[] = {17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 18, 20, 22, 24, 26, 28,
                30, 32, 34, 36, 38, 40, 42, 44, 69};
    for (int i = 0; i < (int)(sizeof(Ls) / sizeof(Ls[0])); i++) {
        int L = Ls[i];
        build_orbit(L);
        printf("L=%2d  p=%ld  N=%ld  p mod stuff below\n", L, pL, Ncp);
        for (int s = 2; s <= 15 && s <= L - 2; s++) analyse(L, L - s);
        free(cps); free(tks); fflush(stdout);
    }
    return 0;
}
