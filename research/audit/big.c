#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline void Astep(int *x, int L) {
    for (int i = 0; i < L - 1; i++) {
        int a = x[i];
        if (!a) continue;
        if (a == 1) { x[i] = 1 + x[i + 1]; x[i + 1] = 0; }
        else { x[i] = a - 1; x[i + 1]++; }
    }
}

static int good(int u, int v) {
    if (u == 0) return v % 2 == 0;
    return (u % 2 == 1) && (v % 2 == 1);
}

/* largest even l such that all aligned pairs among first l coords are good */
static int ell(const int *x, int L) {
    int l = 0;
    while (l + 1 < L && good(x[l], x[l + 1])) l += 2;
    return l;
}

static int zdepth(const int *b, int m) {
    for (int k = 0; k < m; k++) if (b[k] == 0) return k;
    return m;
}

/* ---- full canonical orbit of L=69 ---- */
static void full69(void) {
    const int L = 69;
    int x[69], x0[69];
    for (int i = 0; i < L; i++) x[i] = x0[i] = 1;

    long p = 0, ncp = 0;
    static int gaps[300000];
    static int b0[300000], zb[300000];
    long last = 0;
    for (;;) {
        Astep(x, L); p++;
        if (x[0] == 1) {
            gaps[ncp] = (int)(p - last); last = p;
            b0[ncp] = x[1];
            zb[ncp] = zdepth(x + 1, L - 1);
            ncp++;
        }
        if (!memcmp(x, x0, sizeof x)) break;
    }
    printf("L=69  p_69 = %ld   B_69 = p-68 = %ld\n", p, p - 68);
    printf("      checkpoint period = %ld   Q_69 = %ld\n", ncp, ncp - 1);

    /* gap word palindromic excluding the final R->E transition */
    long Q = ncp;
    int pal_full = 1, pal_drop = 1;
    for (long j = 0; j < Q; j++) if (gaps[j] != gaps[Q - 1 - j]) { pal_full = 0; break; }
    for (long j = 0; j < Q - 1; j++) if (gaps[j] != gaps[Q - 2 - j]) { pal_drop = 0; break; }
    printf("      gap word palindromic (all %ld): %d ; excluding last: %d\n", Q, pal_full, pal_drop);

    /* (b_{Q-j})_0 = z(b_j)  with b_j = H^j(E_68), i.e. gaps/b0/zb indexed from j=1 */
    /* our arrays hold b_1..b_Q where b_Q = b_0 = E_68 ; build index j=0..Q-1 */
    static int B0[300000], ZB[300000];
    B0[0] = b0[Q - 1]; ZB[0] = zb[Q - 1];              /* b_0 = E_68 */
    for (long j = 1; j < Q; j++) { B0[j] = b0[j - 1]; ZB[j] = zb[j - 1]; }
    long QQ = Q - 1; /* Q_69 */
    int idok = 1; long bad = -1;
    for (long j = 0; j <= QQ; j++) {
        long k = ((QQ - j) % Q + Q) % Q;
        if (B0[k] != ZB[j]) { idok = 0; bad = j; break; }
    }
    printf("      (b_{Q-j})_0 = z(b_j) for all j: %d (first failure j=%ld)\n", idok, bad);
}

/* ---- 16-clock outer map for a given L ---- */
static void outer(int L, int maxouter, int verbose) {
    int r = L - 16;
    int *x = malloc(L * sizeof(int));
    for (int i = 0; i < L; i++) x[i] = 1;
    int *y0 = malloc(r * sizeof(int));
    memcpy(y0, x + 16, r * sizeof(int));

    int *w = malloc(r * sizeof(int));
    int *carr = malloc((size_t)maxouter * sizeof(int));
    int *larr = malloc((size_t)maxouter * sizeof(int));

    int cp = 0;                 /* checkpoints since last outer step */
    long steps = 0;
    int nouter = 0, firstescape = -1;
    long period = -1;
    int prefix_ok = 1;

    /* record outer state j=0 (the start) */
    for (;;) {
        /* record current outer state */
        memcpy(w, x + 16, r * sizeof(int));
        for (int t = 0; t < 435; t++) Astep(w, r);
        carr[nouter] = w[0];
        larr[nouter] = ell(x, L);
        for (int i = 0; i < 16; i++) if (x[i] != 1) prefix_ok = 0;
        if (larr[nouter] < 56 && firstescape < 0) firstescape = nouter;
        nouter++;
        if (nouter >= maxouter) break;

        /* advance 112 checkpoints */
        cp = 0;
        while (cp < 112) {
            Astep(x, L); steps++;
            if (x[0] == 1) cp++;
        }
        if (!memcmp(x + 16, y0, r * sizeof(int))) { period = nouter; break; }
    }

    if (verbose) {
        printf("L=%d r=%d : outer period = %ld, prefix E_16 always restored: %d\n",
               L, r, period, prefix_ok);
        if (period > 0) {
            long sumc = 0;
            int cnt[200]; memset(cnt, 0, sizeof cnt);
            for (long j = 0; j < period; j++) { sumc += carr[j]; if (carr[j] < 200) cnt[carr[j]]++; }
            printf("  sum c_j = %ld  -> p = %ld*466 + %ld = %ld   B = %ld\n",
                   sumc, period, sumc, period * 466 + sumc, period * 466 + sumc - (L - 1));
            printf("  c distribution:");
            for (int v = 0; v < 200; v++) if (cnt[v]) printf(" %d:%d", v, cnt[v]);
            printf("\n");
            int lc[200]; memset(lc, 0, sizeof lc);
            int minl = 999;
            for (long j = 0; j < period; j++) { lc[larr[j]]++; if (larr[j] < minl) minl = larr[j]; }
            printf("  min ell = %d ; ell distribution:", minl);
            for (int v = 0; v < 200; v++) if (lc[v]) printf(" %d:%d", v, lc[v]);
            printf("\n");
            /* c word palindromic excluding last */
            int pal = 1;
            for (long j = 0; j < period - 1; j++) if (carr[j] != carr[period - 2 - j]) { pal = 0; break; }
            int palall = 1;
            for (long j = 0; j < period; j++) if (carr[j] != carr[period - 1 - j]) { palall = 0; break; }
            printf("  c word palindromic: all=%d  excluding last=%d\n", palall, pal);
        }
    } else {
        printf("L=%2d  first outer step with ell<56: %s%d   (outer period %ld, examined %d)\n",
               L, firstescape < 0 ? "none up to " : "", firstescape < 0 ? nouter : firstescape,
               period, nouter);
    }
    free(x); free(y0); free(w); free(carr); free(larr);
}

int main(int argc, char **argv) {
    if (argc > 1 && !strcmp(argv[1], "full69")) full69();
    else if (argc > 1 && !strcmp(argv[1], "outer69")) outer(69, 5000, 1);
    else if (argc > 1 && !strcmp(argv[1], "neighbors")) {
        for (int L = 65; L <= 83; L += 2) outer(L, 2000, 0);
    }
    return 0;
}
