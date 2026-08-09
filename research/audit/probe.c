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
static int good(int u, int v) { return u == 0 ? v % 2 == 0 : (u % 2 == 1 && v % 2 == 1); }
static int ell(const int *x, int L) { int l = 0; while (l + 1 < L && good(x[l], x[l + 1])) l += 2; return l; }

/* Is ell>=56 a property of the 16-clock SECTION only, or of every checkpoint? */
static void offsection(void) {
    const int L = 69;
    int x[69], x0[69];
    for (int i = 0; i < L; i++) x[i] = x0[i] = 1;
    int minall = 999, minsec = 999;
    long cnt[70]; memset(cnt, 0, sizeof cnt);
    long ncp = 0;
    for (;;) {
        Astep(x, L);
        if (x[0] == 1) {
            int e = ell(x, L);
            if (e < minall) minall = e;
            cnt[e]++;
            if (ncp % 112 == 111) { if (e < minsec) minsec = e; }
            ncp++;
        }
        if (!memcmp(x, x0, sizeof x)) break;
    }
    printf("L=69: over ALL %ld checkpoints  min ell = %d   (over the 16-clock section: %d)\n",
           ncp, minall, minsec);
    printf("  ell histogram over all checkpoints:");
    for (int e = 0; e <= 68; e++) if (cnt[e]) printf(" %d:%ld", e, cnt[e]);
    printf("\n");
    long below56 = 0; for (int e = 0; e < 56; e++) below56 += cnt[e];
    printf("  checkpoints with ell < 56: %ld  (%.2f%%)\n", below56, 100.0 * below56 / ncp);
}

/* B_L for small L, to check the odd/even growth framing */
static void spectrum(int LMAX) {
    for (int L = 2; L <= LMAX; L++) {
        int *x = malloc(L * sizeof(int)), *x0 = malloc(L * sizeof(int));
        for (int i = 0; i < L; i++) x[i] = x0[i] = 1;
        long p = 0;
        for (;;) { Astep(x, L); p++; if (!memcmp(x, x0, L * sizeof(int))) break; }
        printf("L=%2d  p=%-12ld B=%-12ld %s\n", L, p, p - (L - 1), L % 2 ? "odd" : "even");
        fflush(stdout);
        free(x); free(x0);
    }
}

int main(int argc, char **argv) {
    if (argc > 1 && !strcmp(argv[1], "spectrum")) spectrum(atoi(argv[2]));
    else offsection();
    return 0;
}
