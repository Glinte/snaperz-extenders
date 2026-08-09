/* Round-14 audit, part B: the winching reformulation and the algorithmic eliminations.
     - y-coordinates, L_j gate, A_L = product of L_j
     - L_j = W_j . S_j, order of the ordinary winching sweep O_L
     - defect density, rotating-frame states, branch vectors
     - linear complexity of s_t = (A^t E)_0 mod 2, via LC = n - deg gcd(x^n+1, s(x))  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXL 64

static inline void Astep(int *x, int L) {
    for (int i = 0; i < L - 1; i++) {
        int a = x[i];
        if (!a) continue;
        if (a == 1) { x[i] = 1 + x[i + 1]; x[i + 1] = 0; }
        else { x[i] = a - 1; x[i + 1]++; }
    }
}

/* ---- y coordinates: y_j = L - p_{L-1-j} + j,  p_i = sum_{r<=i} x_r, p_{-1}=0 ---- */
static void to_y(const int *x, int L, int *y) {
    int p[MAXL + 1];
    p[0] = 0;
    for (int i = 0; i < L; i++) p[i + 1] = p[i] + x[i];   /* p[i+1] = p_i */
    for (int j = 0; j <= L; j++) y[j] = L - p[L - j] + j;
}
static void from_y(const int *y, int L, int *x) {
    int p[MAXL + 1];
    for (int j = 0; j <= L; j++) p[L - j] = L - y[j] + j; /* p[i+1] = p_i */
    for (int i = 0; i < L; i++) x[i] = p[i + 1] - p[i];
}

static inline int hj(const int *y, int j) {
    int a = 2 * j, b = y[j + 1] - 1;
    return a < b ? a : b;
}
/* lazy winching gate */
static inline void Lgate(int *y, int j) {
    int d = y[j + 1] - y[j];
    if (d == 1) return;
    if (d == 2) y[j] = y[j - 1] + 1;
    else y[j] = y[j] + 1;
}
/* ordinary bounded winching */
static inline void Wgate(int *y, int j) {
    int h = hj(y, j);
    if (y[j] < h) y[j] = y[j] + 1;
    else y[j] = y[j - 1] + 1;
}
/* top swap of the two highest locally legal values h-1, h */
static inline void Sgate(int *y, int j) {
    int h = hj(y, j);
    if (h - 1 <= y[j - 1]) return;            /* h-1 not legal: act trivially */
    if (y[j] == h) y[j] = h - 1;
    else if (y[j] == h - 1) y[j] = h;
}

static void Lsweep(int *y, int L) { for (int j = L - 1; j >= 1; j--) Lgate(y, j); }
static void Osweep(int *y, int L) { for (int j = L - 1; j >= 1; j--) Wgate(y, j); }

/* ---------- exhaustive small-L structural checks ---------- */
static int Lcur, cnt_states, fail_bij, fail_sweep, fail_fact, fail_sweep_rev;
static int xbuf[MAXL];

static void visit(int L, int *x) {
    int y[MAXL + 1], y2[MAXL + 1], xr[MAXL];
    cnt_states++;
    to_y(x, L, y);
    from_y(y, L, xr);
    if (memcmp(x, xr, L * sizeof(int))) fail_bij++;
    /* sweep in y (j = L-1 .. 1) vs A_L */
    int xa[MAXL];
    memcpy(xa, x, L * sizeof(int));
    Astep(xa, L);
    memcpy(y2, y, (L + 1) * sizeof(int));
    Lsweep(y2, L);
    from_y(y2, L, xr);
    if (memcmp(xa, xr, L * sizeof(int))) fail_sweep++;
    /* the opposite order, j = 1 .. L-1 */
    memcpy(y2, y, (L + 1) * sizeof(int));
    for (int j = 1; j <= L - 1; j++) Lgate(y2, j);
    from_y(y2, L, xr);
    if (memcmp(xa, xr, L * sizeof(int))) fail_sweep_rev++;
    /* L_j = W_j . S_j pointwise, on every state and every j */
    for (int j = 1; j <= L - 1; j++) {
        int ya[MAXL + 1], yb[MAXL + 1];
        memcpy(ya, y, (L + 1) * sizeof(int));
        memcpy(yb, y, (L + 1) * sizeof(int));
        Lgate(ya, j);
        Sgate(yb, j); Wgate(yb, j);
        if (ya[j] != yb[j]) fail_fact++;
    }
}

static void enumerate(int L, int pos, int s) {
    if (pos == L) { if (s == L) visit(L, xbuf); return; }
    int lo = pos + 1 - s; if (lo < 0) lo = 0;
    for (int v = lo; v <= L - s; v++) { xbuf[pos] = v; enumerate(L, pos + 1, s + v); }
    xbuf[pos] = 0;
}

/* ---------- order of the ordinary winching sweep ---------- */
static int order_fail;
static void order_visit(int L, int *x) {
    int y[MAXL + 1], y0[MAXL + 1];
    to_y(x, L, y); memcpy(y0, y, (L + 1) * sizeof(int));
    for (int t = 0; t < 2 * L; t++) Osweep(y, L);
    if (memcmp(y, y0, (L + 1) * sizeof(int))) order_fail++;
}
static void enumerate_ord(int L, int pos, int s) {
    if (pos == L) { if (s == L) order_visit(L, xbuf); return; }
    int lo = pos + 1 - s; if (lo < 0) lo = 0;
    for (int v = lo; v <= L - s; v++) { xbuf[pos] = v; enumerate_ord(L, pos + 1, s + v); }
    xbuf[pos] = 0;
}

/* ---------- 64-bit hashing + distinct count ---------- */
static int cmp64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : (x > y);
}
static long distinct(uint64_t *h, long n) {
    qsort(h, n, sizeof(uint64_t), cmp64);
    long d = n ? 1 : 0;
    for (long i = 1; i < n; i++) if (h[i] != h[i - 1]) d++;
    return d;
}
static inline uint64_t fnv(const void *p, int nbytes) {
    const unsigned char *q = p; uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < nbytes; i++) { h ^= q[i]; h *= 1099511628211ULL; }
    return h;
}

/* ---------- polynomial gcd over F2 for the linear complexity ---------- */
static int deg_of(const uint64_t *a, int nw) {
    for (int w = nw - 1; w >= 0; w--) if (a[w]) { int b = 63; while (!((a[w] >> b) & 1)) b--; return w * 64 + b; }
    return -1;
}
static void xor_shift(uint64_t *a, const uint64_t *b, int nw, int sh) {
    int ws = sh / 64, bs = sh % 64;
    for (int w = nw - 1; w >= ws; w--) {
        uint64_t v = b[w - ws] << bs;
        if (bs && w - ws >= 1) v |= b[w - ws - 1] >> (64 - bs);
        a[w] ^= v;
    }
}
static long linear_complexity(const unsigned char *s, long n) {
    int nw = (int)(n / 64 + 2);
    uint64_t *A = calloc(nw, 8), *B = calloc(nw, 8);
    A[n / 64] |= 1ULL << (n % 64); A[0] ^= 1ULL;                  /* x^n + 1 */
    for (long i = 0; i < n; i++) if (s[i]) B[i / 64] |= 1ULL << (i % 64);
    int da = deg_of(A, nw), db = deg_of(B, nw);
    while (db >= 0) {
        while (da >= db) { xor_shift(A, B, nw, da - db); da = deg_of(A, nw); }
        uint64_t *t = A; A = B; B = t; int td = da; da = db; db = td;
    }
    long g = da < 0 ? 0 : da;
    free(A); free(B);
    return n - g;
}

/* ---------- canonical orbit measurements ---------- */
static void orbit_stats(int L, int do_lc) {
    int x[MAXL], x0[MAXL];
    for (int i = 0; i < L; i++) x[i] = x0[i] = 1;
    /* first pass: period */
    long p = 0;
    for (;;) { Astep(x, L); p++; if (!memcmp(x, x0, L * sizeof(int))) break; }

    uint64_t *hr = malloc(p * 8), *hb = malloc(p * 8);
    unsigned char *bits = malloc(p);
    long defect = 0, gates = 0, sweeps_with = 0;
    long diffcoord = 0;

    for (int i = 0; i < L; i++) x[i] = 1;
    int yrot[MAXL + 1];
    for (long t = 0; t < p; t++) {
        /* branch vector of this sweep + defect count, from the y-picture */
        int y[MAXL + 1], yl[MAXL + 1], yw[MAXL + 1];
        to_y(x, L, y);
        unsigned char br[MAXL];
        for (int i = 0; i < L - 1; i++) br[i] = x[i] == 0 ? 0 : (x[i] == 1 ? 1 : 2);
        hb[t] = fnv(br, L - 1);

        memcpy(yl, y, (L + 1) * sizeof(int));
        int local = 0;
        for (int j = L - 1; j >= 1; j--) {
            memcpy(yw, yl, (L + 1) * sizeof(int));
            int before = yl[j];
            Lgate(yl, j);
            Wgate(yw, j);
            if (yl[j] != yw[j]) { defect++; local++; }
            (void)before;
            gates++;
        }
        if (local) sweeps_with++;

        /* rotating frame r_t = O^{-t} A^t (E) */
        to_y(x, L, yrot);
        long back = (2L * L - (t % (2L * L))) % (2L * L);
        for (long u = 0; u < back; u++) Osweep(yrot, L);
        int xr[MAXL];
        from_y(yrot, L, xr);
        hr[t] = fnv(xr, L * sizeof(int));
        for (int i = 0; i < L; i++) if (xr[i] != 1) diffcoord++;

        bits[t] = (unsigned char)(x[0] & 1);
        Astep(x, L);
    }
    long dr = distinct(hr, p), dbv = distinct(hb, p);
    printf("L=%2d  p_L=%-8ld  defective gates %.2f%%  defects/sweep %.3f  "
           "rot-frame distinct %ld  branch vectors distinct %ld (%.1f%% unique-ish)  "
           "mean coords != 1 in rot frame %.2f\n",
           L, p, 100.0 * defect / gates, (double)defect / p, dr, dbv,
           100.0 * dbv / p, (double)diffcoord / p);
    if (do_lc) {
        long lc = linear_complexity(bits, p);
        printf("      linear complexity of (A^t E)_0 mod 2 over F2: %ld  (p_L = %ld)\n", lc, p);
    }
    free(hr); free(hb); free(bits);
}

int main(int argc, char **argv) {
    printf("--- exhaustive structural checks ---\n");
    for (int L = 5; L <= 10; L++) {
        cnt_states = fail_bij = fail_sweep = fail_fact = fail_sweep_rev = 0;
        enumerate(L, 0, 0);
        order_fail = 0;
        enumerate_ord(L, 0, 0);
        printf("L=%2d |C_L|=%6d  y-bijection %s  A_L = L_{L-1}..L_1 (j desc) %s  (j asc) %s  "
               "L_j = W_j.S_j %s (%d fails)  O_L^{2L}=id %s\n",
               L, cnt_states, fail_bij ? "FAIL" : "OK", fail_sweep ? "FAIL" : "OK",
               fail_sweep_rev ? "FAIL" : "OK", fail_fact ? "FAIL" : "OK", fail_fact,
               order_fail ? "FAIL" : "OK");
    }
    printf("\n--- canonical orbit measurements ---\n");
    int Ls[] = {25, 29, 31, 35, 36};
    for (unsigned i = 0; i < sizeof(Ls) / sizeof(Ls[0]); i++)
        orbit_stats(Ls[i], argc > 1 ? 1 : 0);
    return 0;
}
