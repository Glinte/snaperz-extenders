/* Round-27b: the exceptional-ECO edge set at L=69, and what pins the 56|13
 * excursion boundaries to it.
 *
 * r27_reflect showed the clean 56|13 checkpoint section is *totally* disjoint
 * from its own reflection (0 of 268,128), so the 84-excursion palindrome cannot
 * come from that section being sigma-invariant.  But the excursion intervals
 * still land on their partners to within a small offset.
 *
 * The reversor does provably act on the exceptional-ECO edge set: e_s = e_{-s}.
 * So the natural bridge is that each excursion boundary sits at a bounded,
 * reflection-covariant distance from an exceptional edge.  This measures:
 *
 *   1. e_t as an indicator, and its symmetry under both the state centre
 *      sigma(t) = C - t and the edge centre C - 1 - t (the correct one for
 *      edges is the question, so both are reported);
 *   2. the distance from each excursion boundary to the nearest exceptional
 *      edge, and whether those distances are themselves reflected.
 *
 * usage: r27b_eco [L] [cut]      (defaults 69 56)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int L, CUT, S;
static uint8_t *cyc;
static int pS;

static void sweep_n(uint8_t *x, int n) {
    for (int i = 0; i < n - 1; i++) {
        uint8_t a = x[i];
        if (!a) continue;
        if (a == 1) { x[i] = 1 + x[i + 1]; x[i + 1] = 0; }
        else { x[i] = a - 1; x[i + 1] += 1; }
    }
}

static int last_pos(const uint8_t *x, int n) {
    int m = n - 1;
    while (m >= 0 && !x[m]) m--;
    return m;
}
static int hgt(const uint8_t *x, int n) { return n - last_pos(x, n); }

/* unlam: delete the rightmost peak of a length-n state -> (t of length n-1, r) */
static void unlam(const uint8_t *x, int n, uint8_t *t, int *r) {
    int m = last_pos(x, n);
    if (x[m] == 1) {
        for (int i = 0; i < m; i++) t[i] = x[i];
        for (int i = m; i < n - 1; i++) t[i] = 0;
        *r = n - 1 - m;
    } else {
        for (int i = 0; i < m; i++) t[i] = x[i];
        t[m] = x[m] - 1;
        for (int i = m + 1; i < n - 1; i++) t[i] = 0;
        *r = hgt(t, n - 1);
    }
}

/* edge out of x is exceptional iff h(A t) = h(t) - 1 and r = h(t) */
static int is_exceptional(const uint8_t *x, uint8_t *t, uint8_t *at) {
    int r;
    unlam(x, L, t, &r);
    int ht = hgt(t, L - 1);
    memcpy(at, t, L - 1);
    sweep_n(at, L - 1);
    return (hgt(at, L - 1) == ht - 1) && (r == ht);
}

static void build_cycle(void) {
    uint8_t *x = calloc(S, 1), *x0 = calloc(S, 1);
    for (int i = 0; i < S; i++) x0[i] = 1;
    memcpy(x, x0, S);
    cyc = malloc((size_t)(1 << 16) * S);
    pS = 0;
    do { memcpy(cyc + (size_t)pS * S, x, S); pS++; sweep_n(x, S); }
    while (memcmp(x, x0, S));
    free(x); free(x0);
}
static int phase_of(const uint8_t *tail) {
    for (int k = 0; k < pS; k++)
        if (!memcmp(cyc + (size_t)k * S, tail, S)) return k;
    return -1;
}

int main(int argc, char **argv) {
    L = argc > 1 ? atoi(argv[1]) : 69;
    CUT = argc > 2 ? atoi(argv[2]) : 56;
    S = L - CUT;
    build_cycle();

    uint8_t *x = calloc(L, 1), *x0 = calloc(L, 1);
    uint8_t *t = calloc(L, 1), *at = calloc(L, 1);
    for (int i = 0; i < L; i++) x0[i] = 1;

    memcpy(x, x0, L);
    long p = 0;
    do { sweep_n(x, L); p++; } while (memcmp(x, x0, L));

    uint8_t *exc = calloc(p, 1), *is_clean = calloc(p, 1);
    int32_t *dat = malloc(p * sizeof(int32_t));
    for (long i = 0; i < p; i++) dat[i] = -1;

    memcpy(x, x0, L);
    long nexc = 0, nclean = 0;
    for (long s = 0; s < p; s++) {
        if (is_exceptional(x, t, at)) { exc[s] = 1; nexc++; }
        if (x[0] == 1) {
            int sm = 0, clean = 1;
            for (int i = 0; i < CUT; i++) { sm += x[i]; if (sm < i + 1) clean = 0; }
            if (sm != CUT) clean = 0;
            if (clean) {
                int ph = phase_of(x + CUT);
                if (ph >= 0) {
                    is_clean[s] = 1; nclean++;
                    dat[s] = (int32_t)((((long)ph - s) % pS + pS) % pS);
                }
            }
        }
        sweep_n(x, L);
    }

    long C = ((p - 2 * (long)L + 2) % p + p) % p;
    printf("L=%d cut=%d|%d  p_L=%ld  p_%d=%d  centre C=%ld\n", L, CUT, S, p, S, pS, C);
    printf("exceptional edges: %ld  (Xi_period at L=69 is 29,265)\n", nexc);

    long sym_state = 0, sym_edge = 0;
    for (long s = 0; s < p; s++) {
        if (!exc[s]) continue;
        long u1 = ((C - s) % p + p) % p;
        long u2 = ((C - 1 - s) % p + p) % p;
        if (exc[u1]) sym_state++;
        if (exc[u2]) sym_edge++;
    }
    printf("  e_s = e_{C-s}   : %ld/%ld\n", sym_state, nexc);
    printf("  e_s = e_{C-1-s} : %ld/%ld\n", sym_edge, nexc);

    /* excursion boundaries */
    long *ct = malloc(nclean * sizeof(long)), m = 0;
    for (long s = 0; s < p; s++) if (is_clean[s]) ct[m++] = s;
    long nex = 0;
    for (long i = 1; i < m; i++) if (dat[ct[i]] != dat[ct[i - 1]]) nex++;
    long *b0 = malloc(nex * sizeof(long)), *b1 = malloc(nex * sizeof(long));
    long e = 0;
    for (long i = 1; i < m; i++)
        if (dat[ct[i]] != dat[ct[i - 1]]) { b0[e] = ct[i - 1]; b1[e] = ct[i]; e++; }

    /* distance from each boundary forward to the next exceptional edge */
    long *fwd = malloc(p * sizeof(long));
    long d = 0;
    for (long s = p - 1; s >= 0; s--) { fwd[s] = exc[s] ? 0 : d + 1; d = fwd[s]; }
    for (int pass = 0; pass < 2; pass++)
        for (long s = p - 1; s >= 0; s--) {
            long nx = (s + 1) % p;
            if (!exc[s] && fwd[nx] + 1 < fwd[s]) fwd[s] = fwd[nx] + 1;
        }

    printf("\n  distance from excursion boundary to next exceptional edge:\n");
    long mx0 = 0, mx1 = 0, refl_ok = 0;
    for (long i = 0; i < nex; i++) {
        if (fwd[b0[i]] > mx0) mx0 = fwd[b0[i]];
        if (fwd[b1[i]] > mx1) mx1 = fwd[b1[i]];
        long k = nex - 1 - i;
        if (fwd[b0[i]] == fwd[b1[k]] || fwd[b1[i]] == fwd[b0[k]]) refl_ok++;
    }
    printf("    max over b0: %ld   max over b1: %ld   (excursions are 358..89,214 long)\n",
           mx0, mx1);
    printf("    boundaries whose distance matches the reflected partner: %ld/%ld\n",
           refl_ok, nex);
    printf("    first 8 (b0,dist) (b1,dist): ");
    for (long i = 0; i < 8 && i < nex; i++)
        printf("(%ld,%ld)(%ld,%ld) ", b0[i], fwd[b0[i]], b1[i], fwd[b1[i]]);
    printf("\n");

    /* The interval is NOT exactly reflected (76/83 lengths differ), but the
     * exceptional-edge set is exactly sigma-invariant.  So count exceptional
     * edges per excursion and ask whether *that* is palindromic. */
    long *pre = malloc((p + 1) * sizeof(long));
    pre[0] = 0;
    for (long s = 0; s < p; s++) pre[s + 1] = pre[s] + exc[s];

    long *cnt = malloc(nex * sizeof(long));
    for (long i = 0; i < nex; i++) cnt[i] = pre[b1[i]] - pre[b0[i]];

    long pal = 0;
    for (long i = 0; i < nex; i++) if (cnt[i] == cnt[nex - 1 - i]) pal++;
    printf("\n  exceptional edges strictly inside each excursion:\n");
    printf("    palindromic (cnt[e] == cnt[nex-1-e]): %ld/%ld\n", pal, nex);
    printf("    counts: ");
    for (long i = 0; i < nex && i < 30; i++) printf("%ld ", cnt[i]);
    printf("%s\n", nex > 30 ? "..." : "");

    /* does the jump depend only on that count? */
    printf("    (count -> jump) map, checking it is single-valued:\n");
    long uc[64], uj[64], nu = 0;
    int multi = 0;
    for (long i = 0; i < nex; i++) {
        long jump = 0;
        /* recover the jump from dat[] at the two boundaries */
        jump = (((long)dat[b1[i]] - dat[b0[i]]) % pS + pS) % pS;
        int f = -1;
        for (long j = 0; j < nu; j++) if (uc[j] == cnt[i]) f = (int)j;
        if (f < 0 && nu < 64) { f = (int)nu++; uc[f] = cnt[i]; uj[f] = jump; }
        else if (f >= 0 && uj[f] != jump) multi = 1;
    }
    for (long j = 0; j < nu && j < 20; j++)
        printf("      count=%-8ld jump=%ld\n", uc[j], uj[j]);
    printf("    single-valued: %s\n", multi ? "NO" : "YES");
    return 0;
}
