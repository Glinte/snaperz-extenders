/* Round-14 audit, part C: pin down the L_j = W_j.S_j factorization and
   redo the branch-vector count with the branch read *during* the sweep. */
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
/* branch vector read at the moment each gate fires */
static inline void Astep_br(int *x, int L, unsigned char *br) {
    for (int i = 0; i < L - 1; i++) {
        int a = x[i];
        br[i] = a == 0 ? 0 : (a == 1 ? 1 : 2);
        if (!a) continue;
        if (a == 1) { x[i] = 1 + x[i + 1]; x[i + 1] = 0; }
        else { x[i] = a - 1; x[i + 1]++; }
    }
}

static void to_y(const int *x, int L, int *y) {
    int p[MAXL + 1];
    p[0] = 0;
    for (int i = 0; i < L; i++) p[i + 1] = p[i] + x[i];
    for (int j = 0; j <= L; j++) y[j] = L - p[L - j] + j;
}
static inline int hj(const int *y, int j) {
    int a = 2 * j, b = y[j + 1] - 1;
    return a < b ? a : b;
}
static inline int Lval(const int *y, int j) {
    int d = y[j + 1] - y[j];
    if (d == 1) return y[j];
    if (d == 2) return y[j - 1] + 1;
    return y[j] + 1;
}
static inline int Wval(const int *y, int j) {
    int h = hj(y, j);
    return y[j] < h ? y[j] + 1 : y[j - 1] + 1;
}
/* S then W, with the "both top values legal" guard */
static inline int WSval(const int *y, int j) {
    int h = hj(y, j), t[MAXL + 1];
    memcpy(t, y, (j + 2) * sizeof(int));
    if (h - 1 > y[j - 1]) {
        if (t[j] == h) t[j] = h - 1;
        else if (t[j] == h - 1) t[j] = h;
    }
    return Wval(t, j);
}
/* W then S */
static inline int SWval(const int *y, int j) {
    int t[MAXL + 1];
    memcpy(t, y, (j + 2) * sizeof(int));
    t[j] = Wval(y, j);
    int h = hj(t, j);
    if (h - 1 > t[j - 1]) {
        if (t[j] == h) t[j] = h - 1;
        else if (t[j] == h - 1) t[j] = h;
    }
    return t[j];
}

static int xbuf[MAXL];
static long n_gate, n_defect, n_ws_fail, n_sw_fail, n_ceiling_fail, n_other;
static int shown;

static void visit(int L, int *x) {
    int y[MAXL + 1];
    to_y(x, L, y);
    for (int j = 1; j <= L - 1; j++) {
        n_gate++;
        int lv = Lval(y, j), wv = Wval(y, j), ws = WSval(y, j), sw = SWval(y, j);
        if (lv != wv) n_defect++;
        if (lv != ws) {
            n_ws_fail++;
            int atceil = (y[j] == 2 * j);
            if (atceil) n_ceiling_fail++; else n_other++;
            if (shown < 6) {
                printf("   fail L=%d j=%d  y[j-1..j+1]=(%d,%d,%d)  h=%d  2j=%d  "
                       "L_j=%d  (W.S)_j=%d  W_j=%d  at-ceiling=%d\n",
                       L, j, y[j - 1], y[j], y[j + 1], hj(y, j), 2 * j, lv, ws, wv, atceil);
                shown++;
            }
        }
        if (lv != sw) n_sw_fail++;
    }
}
static void enumerate(int L, int pos, int s) {
    if (pos == L) { if (s == L) visit(L, xbuf); return; }
    int lo = pos + 1 - s; if (lo < 0) lo = 0;
    for (int v = lo; v <= L - s; v++) { xbuf[pos] = v; enumerate(L, pos + 1, s + v); }
    xbuf[pos] = 0;
}

static int cmp64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : (x > y);
}
static inline uint64_t fnv(const void *p, int n) {
    const unsigned char *q = p; uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < n; i++) { h ^= q[i]; h *= 1099511628211ULL; }
    return h;
}

static void branch_count(int L) {
    int x[MAXL], x0[MAXL];
    for (int i = 0; i < L; i++) x[i] = x0[i] = 1;
    long p = 0;
    for (;;) { Astep(x, L); p++; if (!memcmp(x, x0, L * sizeof(int))) break; }
    uint64_t *h = malloc(p * 8);
    unsigned char br[MAXL];
    for (int i = 0; i < L; i++) x[i] = 1;
    long consec = 0;
    uint64_t prev = 0;
    for (long t = 0; t < p; t++) {
        Astep_br(x, L, br);
        h[t] = fnv(br, L - 1);
        if (t && h[t] == prev) consec++;
        prev = h[t];
    }
    qsort(h, p, sizeof(uint64_t), cmp64);
    long d = 1;
    for (long i = 1; i < p; i++) if (h[i] != h[i - 1]) d++;
    printf("L=%2d  p_L=%-8ld  distinct branch vectors (read during the sweep) = %ld  "
           "(%.1f%%)  consecutive repeats = %ld\n", L, p, d, 100.0 * d / p, consec);
    free(h);
}

int main(void) {
    printf("--- L_j = W_j . S_j, exhaustive ---\n");
    for (int L = 5; L <= 10; L++) {
        n_gate = n_defect = n_ws_fail = n_sw_fail = n_ceiling_fail = n_other = 0;
        enumerate(L, 0, 0);
        printf("L=%2d  gates %ld  defects(L!=W) %ld (%.2f%%)  L!=(W.S) %ld (%.2f%%)  "
               "L!=(S.W) %ld   of the W.S failures: at y_j=2j %ld, elsewhere %ld\n",
               L, n_gate, n_defect, 100.0 * n_defect / n_gate, n_ws_fail,
               100.0 * n_ws_fail / n_gate, n_sw_fail, n_ceiling_fail, n_other);
    }
    printf("\n--- branch vectors, read during the sweep ---\n");
    int Ls[] = {25, 29, 31, 35};
    for (unsigned i = 0; i < sizeof(Ls) / sizeof(Ls[0]); i++) branch_count(Ls[i]);
    return 0;
}
