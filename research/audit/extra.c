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

/* trace ell around the claimed reflection at outer step 780/781 */
static void trace(int L, int lo, int hi) {
    int r = L - 16, *x = malloc(L * sizeof(int));
    for (int i = 0; i < L; i++) x[i] = 1;
    printf("L=%d ell(x_j) for j=%d..%d:", L, lo, hi);
    for (int j = 0; j <= hi; j++) {
        if (j >= lo) printf(" %d", ell(x, L));
        int cp = 0; while (cp < 112) { Astep(x, L); if (x[0] == 1) cp++; }
    }
    printf("\n"); free(x); (void)r;
}

/* does a "clock" exist at prefix size P?  i.e. is there k<=KMAX with
   H^k(E_{P}, y) having prefix E_P, simultaneously for all small y ? */
static void clock_test(int P, int KMAX) {
    /* full state is (E_P, y) of length P+r, checkpoints are x[0]==1 */
    printf("prefix P=%2d : ", P);
    int found = -1;
    for (int k = 1; k <= KMAX && found < 0; k++) {
        int ok = 1;
        for (int r = 1; r <= 4 && ok; r++) {
            /* test a couple of suffixes: E_r and one other state */
            for (int variant = 0; variant < 2 && ok; variant++) {
                int L = P + r, *x = malloc(L * sizeof(int));
                for (int i = 0; i < L; i++) x[i] = 1;
                if (variant == 1) { for (int t = 0; t < 3; t++) Astep(x + P, r); }
                if (x[0] != 1) { free(x); continue; }
                int cp = 0;
                while (cp < k) { Astep(x, L); if (x[0] == 1) cp++; }
                for (int i = 0; i < P; i++) if (x[i] != 1) ok = 0;
                free(x);
            }
        }
        if (ok) found = k;
    }
    if (found > 0) printf("clock found, k=%d\n", found);
    else printf("no clock with k<=%d\n", KMAX);
}

int main(int argc, char **argv) {
    if (argc > 1 && !strcmp(argv[1], "trace")) {
        trace(69, 775, 790);
        for (int L = 65; L <= 75; L += 2) if (L != 69) trace(L, 775, 790);
    } else {
        for (int P = 2; P <= 64; P *= 2) clock_test(P, 3000);
    }
    return 0;
}
