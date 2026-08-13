/* Round-26 audit: eps_L as a midpoint observable.
 *
 * Claim under test:  eps_L = 1  iff  p_L is even and h(A^{p_L/2} R_L) = 1,
 * hence T_L = (1 + eps_L) p_L - (L - 1).  Reflection pairs every non-fixed
 * phase of the canonical cycle, so the winner bit is decided at the only other
 * fixed point of the reversor E_L -- the midpoint.
 *
 * Walks the canonical cycle from R_L = (L, 0, ..., 0) using the bare gate rule,
 * so this is independent of the Python core.  Compare the output against the
 * community T_L table (see docs/research-notes.md, Data sources).
 *
 * usage: r26_eps Lmin Lmax     ->  L p_L p_even eps T_L   per line
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int L;

static void sweep(unsigned char *x) {
    for (int i = 0; i < L - 1; i++) {
        unsigned char a = x[i];
        if (!a) continue;
        if (a == 1) { x[i] = 1 + x[i + 1]; x[i + 1] = 0; }
        else        { x[i] = a - 1; x[i + 1] += 1; }
    }
}

static int hgt(const unsigned char *x) {
    int m = L - 1;
    while (m >= 0 && !x[m]) m--;
    return L - m;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s Lmin Lmax\n", argv[0]); return 1; }
    int lo = atoi(argv[1]), hi = atoi(argv[2]);

    for (L = lo; L <= hi; L++) {
        unsigned char *R = calloc(L, 1), *x = calloc(L, 1);
        R[0] = (unsigned char)L;
        memcpy(x, R, L);

        long long p = 0;
        do { sweep(x); p++; } while (memcmp(x, R, L));

        int hmid = -1;
        if (p % 2 == 0) {
            memcpy(x, R, L);
            for (long long i = 0; i < p / 2; i++) sweep(x);
            hmid = hgt(x);
        }
        int eps = (p % 2 == 0 && hmid == 1) ? 1 : 0;
        printf("%d %lld %d %d %lld\n", L, p, (int)(p % 2 == 0), eps,
               (long long)(1 + eps) * p - (L - 1));
        free(R); free(x);
    }
    return 0;
}
