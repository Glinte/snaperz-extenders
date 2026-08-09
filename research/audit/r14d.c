/* Round-14 audit, part D:
     - the repaired factorization L_j = W_j . S_j (swap only when the ceiling is slack)
     - linear complexity of s_t = (A^t E)_0 mod 2
     - the L=69 separator word: palindrome decomposition and linear complexity  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXL 80

static inline void Astep(int *x, int L) {
    for (int i = 0; i < L - 1; i++) {
        int a = x[i];
        if (!a) continue;
        if (a == 1) { x[i] = 1 + x[i + 1]; x[i + 1] = 0; }
        else { x[i] = a - 1; x[i + 1]++; }
    }
}
static void to_y(const int *x, int L, int *y) {
    int p[MAXL + 1]; p[0] = 0;
    for (int i = 0; i < L; i++) p[i + 1] = p[i] + x[i];
    for (int j = 0; j <= L; j++) y[j] = L - p[L - j] + j;
}
static inline int hj(const int *y, int j) { int a = 2 * j, b = y[j + 1] - 1; return a < b ? a : b; }
static inline int Lval(const int *y, int j) {
    int d = y[j + 1] - y[j];
    if (d == 1) return y[j];
    if (d == 2) return y[j - 1] + 1;
    return y[j] + 1;
}
static inline int Wval_at(const int *y, int j, int v) {
    int h = hj(y, j);
    return v < h ? v + 1 : y[j - 1] + 1;
}
/* repaired: apply the top swap only when the staircase ceiling is NOT the binding
   constraint, i.e. only when h_j = y_{j+1}-1 (equivalently 2j >= y_{j+1}-1). */
static inline int WSval_fixed(const int *y, int j) {
    int h = hj(y, j), v = y[j];
    int ceiling_slack = (2 * j >= y[j + 1] - 1);
    if (ceiling_slack && h - 1 > y[j - 1]) {
        if (v == h) v = h - 1;
        else if (v == h - 1) v = h;
    }
    return Wval_at(y, j, v);
}

static int xbuf[MAXL];
static long g_gate, g_fail;
static void visit(int L, int *x) {
    int y[MAXL + 1]; to_y(x, L, y);
    for (int j = 1; j <= L - 1; j++) { g_gate++; if (Lval(y, j) != WSval_fixed(y, j)) g_fail++; }
}
static void enumerate(int L, int pos, int s) {
    if (pos == L) { if (s == L) visit(L, xbuf); return; }
    int lo = pos + 1 - s; if (lo < 0) lo = 0;
    for (int v = lo; v <= L - s; v++) { xbuf[pos] = v; enumerate(L, pos + 1, s + v); }
    xbuf[pos] = 0;
}

/* ---- F2 polynomial gcd -> linear complexity of a periodic sequence ---- */
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
    A[n / 64] |= 1ULL << (n % 64); A[0] ^= 1ULL;
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

static void lc_of_orbit(int L) {
    int x[MAXL], x0[MAXL];
    for (int i = 0; i < L; i++) x[i] = x0[i] = 1;
    long p = 0;
    for (;;) { Astep(x, L); p++; if (!memcmp(x, x0, L * sizeof(int))) break; }
    unsigned char *bits = malloc(p);
    for (int i = 0; i < L; i++) x[i] = 1;
    for (long t = 0; t < p; t++) { bits[t] = (unsigned char)(x[0] & 1); Astep(x, L); }
    printf("L=%2d  p_L=%-8ld  LC of (A^t E)_0 mod 2 = %ld\n", L, p, linear_complexity(bits, p));
    free(bits);
}

/* ---- L=69 separator word ---- */
static int sep[3000];
static void sep69(void) {
    const int L = 69;
    int x[69], x0[69];
    for (int i = 0; i < L; i++) x[i] = x0[i] = 1;
    long p = 0, k = 0, last = 0; int n = 0;
    for (;;) {
        Astep(x, L); p++;
        if (x[0] == 1) {
            int g = (int)(p - last); last = p;
            if (k % 112 == 111) sep[n++] = g;
            k++;
        }
        if (!memcmp(x, x0, sizeof x)) break;
    }
    printf("\n[69] separator word length %d\n", n);
    int cnt[128]; memset(cnt, 0, sizeof cnt);
    for (int i = 0; i < n; i++) if (sep[i] < 128) cnt[sep[i]]++;
    printf("     distribution:");
    for (int v = 0; v < 128; v++) if (cnt[v]) printf(" %d^%d", v, cnt[v]);
    printf("   (claim 18^1050 20^718 22^395 24^158 26^58 28^14 69^1)\n");

    int m = n - 1;                       /* first n-1 symbols */
    int pal = 1;
    for (int i = 0; i < m; i++) if (sep[i] != sep[m - 1 - i]) { pal = 0; break; }
    printf("     first %d symbols palindromic: %s\n", m, pal ? "OK" : "FAIL");
    /* claimed decomposition P_935 M_523 P_935 with P, M themselves palindromes */
    if (m == 2 * 935 + 523) {
        int okP = 1, okM = 1;
        for (int i = 0; i < 935; i++) if (sep[i] != sep[934 - i]) { okP = 0; break; }
        for (int i = 0; i < 523; i++) if (sep[935 + i] != sep[935 + 522 - i]) { okM = 0; break; }
        int okrep = !memcmp(sep, sep + 935 + 523, 935 * sizeof(int));
        printf("     P_935 M_523 P_935: P palindrome %s, M palindrome %s, both P blocks equal %s\n",
               okP ? "OK" : "FAIL", okM ? "OK" : "FAIL", okrep ? "OK" : "FAIL");
    } else printf("     length %d does not match 935+523+935=%d\n", m, 2 * 935 + 523);

    unsigned char *b = malloc(n);
    for (int i = 0; i < n; i++) b[i] = sep[i] == 20;
    printf("     LC[s_j=20] = %ld  (claim 2392)\n", linear_complexity(b, n));
    for (int i = 0; i < n; i++) b[i] = sep[i] >= 24;
    printf("     LC[s_j>=24] = %ld  (claim 2390)\n", linear_complexity(b, n));
    for (int i = 0; i < n; i++) b[i] = sep[i] & 1;
    printf("     LC[parity]  = %ld  (claim 2394)\n", linear_complexity(b, n));
    free(b);
}

int main(void) {
    printf("--- repaired factorization: swap only when the ceiling is slack ---\n");
    for (int L = 5; L <= 11; L++) {
        g_gate = g_fail = 0;
        enumerate(L, 0, 0);
        printf("L=%2d  gates %ld  L_j != W_j.S_j : %ld  %s\n", L, g_gate, g_fail, g_fail ? "FAIL" : "OK");
    }
    printf("\n--- linear complexity of the canonical root-parity sequence ---\n");
    int Ls[] = {25, 29, 31, 35};
    for (unsigned i = 0; i < sizeof(Ls) / sizeof(Ls[0]); i++) lc_of_orbit(Ls[i]);
    sep69();
    return 0;
}
