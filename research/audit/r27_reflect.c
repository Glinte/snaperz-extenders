/* Round-27: does the pointed reversor explain the L=69 excursion palindrome?
 *
 * The reversor gives E_L(A^s R_L) = A^{-s} R_L.  Rooted at E_L = A^{L-1} R_L,
 * that is the involution on times
 *
 *     sigma(t) = (p_L - 2L + 2) - t   (mod p_L),
 *
 * the same centre the global height reflection delta_t = -delta_{p_L-2L+1-t}
 * uses.  The 84-excursion palindrome at the clean 56|13 cut is the last L=69
 * conjecture the reversor has not been shown to cover: it proves the canonical
 * height and exceptional-ECO palindromes, but the 56|13 defect classification
 * is a *different* section (a_0 = 1 plus a clean prefix), and nothing so far
 * says that section is sigma-invariant.
 *
 * This measures exactly that, in increasing order of specificity:
 *   1. is the checkpoint-time set sigma-invariant?
 *   2. is the clean-checkpoint subset sigma-invariant?
 *   3. is the set of constant-d run boundaries sigma-invariant?
 *   4. does excursion e map to excursion (nex-1-e), i.e. is the palindrome
 *      literally the reversor acting on the excursion sequence?
 *
 * usage: r27_reflect [L] [cut]        (defaults 69 56)
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
    for (int i = 0; i < L; i++) x0[i] = 1;
    memcpy(x, x0, L);

    /* first pass: p_L */
    long p = 0;
    do { sweep_n(x, L); p++; } while (memcmp(x, x0, L));

    uint8_t *is_cp = calloc(p, 1), *is_clean = calloc(p, 1);
    int32_t *dat = malloc(p * sizeof(int32_t));
    for (long t = 0; t < p; t++) dat[t] = -1;

    memcpy(x, x0, L);
    long ncp = 0, nclean = 0;
    for (long t = 0; t < p; t++) {
        if (x[0] == 1) {
            is_cp[t] = 1; ncp++;
            int s = 0, clean = 1;
            for (int i = 0; i < CUT; i++) { s += x[i]; if (s < i + 1) clean = 0; }
            if (s != CUT) clean = 0;
            if (clean) {
                int ph = phase_of(x + CUT);
                if (ph >= 0) {
                    is_clean[t] = 1; nclean++;
                    dat[t] = (int32_t)((((long)ph - t) % pS + pS) % pS);
                }
            }
        }
        sweep_n(x, L);
    }

    /* sigma(t) = C - t mod p.  Note 2L would be the long literal 2 in C. */
    long C = ((p - 2 * (long)L + 2) % p + p) % p;
    printf("L=%d cut=%d|%d  p_L=%ld  p_%d=%d\n", L, CUT, S, p, S, pS);
    printf("reflection centre: sigma(t) = %ld - t  (mod %ld)\n\n", C, p);

    long cp_sym = 0, clean_sym = 0;
    for (long t = 0; t < p; t++) {
        long u = ((C - t) % p + p) % p;
        if (is_cp[t] && is_cp[u]) cp_sym++;
        if (is_clean[t] && is_clean[u]) clean_sym++;
    }
    printf("1. checkpoints       : %ld total, %ld with sigma-partner also a checkpoint (%.4f%%)\n",
           ncp, cp_sym, 100.0 * cp_sym / ncp);
    printf("2. clean on-cycle    : %ld total, %ld with sigma-partner also clean    (%.4f%%)\n",
           nclean, clean_sym, 100.0 * clean_sym / nclean);

    /* run boundaries over clean on-cycle checkpoints in time order */
    long *ct = malloc(nclean * sizeof(long));
    long m = 0;
    for (long t = 0; t < p; t++) if (is_clean[t]) ct[m++] = t;

    long nex = 0;
    for (long i = 1; i < m; i++) if (dat[ct[i]] != dat[ct[i - 1]]) nex++;
    long *b0 = malloc((nex + 1) * sizeof(long));  /* last cp of run e   */
    long *b1 = malloc((nex + 1) * sizeof(long));  /* first cp of run e+1 */
    int  *jm = malloc((nex + 1) * sizeof(int));
    long e = 0;
    for (long i = 1; i < m; i++)
        if (dat[ct[i]] != dat[ct[i - 1]]) {
            b0[e] = ct[i - 1]; b1[e] = ct[i];
            jm[e] = (int)(((dat[ct[i]] - dat[ct[i - 1]]) % pS + pS) % pS);
            e++;
        }
    printf("3. excursions        : %ld\n", nex);

    /* is the boundary set sigma-invariant, and does e <-> nex-1-e? */
    long hit_b0 = 0, hit_b1 = 0, cross = 0, jrev = 0;
    for (long i = 0; i < nex; i++) {
        long s0 = ((C - b0[i]) % p + p) % p;
        long s1 = ((C - b1[i]) % p + p) % p;
        for (long j = 0; j < nex; j++) {
            if (b0[j] == s0 || b1[j] == s0) { hit_b0++; break; }
        }
        for (long j = 0; j < nex; j++) {
            if (b0[j] == s1 || b1[j] == s1) { hit_b1++; break; }
        }
        long k = nex - 1 - i;
        if (s1 >= b0[k] && s1 <= b1[k]) cross++;
        if (jm[i] == jm[k]) jrev++;
    }
    printf("   boundary sigma-images landing on a boundary: b0 %ld/%ld, b1 %ld/%ld\n",
           hit_b0, nex, hit_b1, nex);
    printf("   sigma(b1_e) inside excursion (nex-1-e)      : %ld/%ld\n", cross, nex);
    printf("4. jump word palindromic (jm[e]==jm[nex-1-e])  : %ld/%ld\n", jrev, nex);

    /* where do the sigma-images actually land relative to the boundaries? */
    /* offsets sigma(b1_e) - b0_k and sigma(b0_e) - b1_k, over all e */
    printf("\n   offsets of the reflected excursion against its partner:\n");
    long off_lo[8], off_hi[8], nlo = 0, nhi = 0, exact = 0, contained = 0;
    for (long i = 0; i < nex; i++) {
        long k = nex - 1 - i;
        long s0 = ((C - b0[i]) % p + p) % p;
        long s1 = ((C - b1[i]) % p + p) % p;
        long dlo = s1 - b0[k], dhi = s0 - b1[k];
        int f = 0;
        for (long j = 0; j < nlo; j++) if (off_lo[j] == dlo) f = 1;
        if (!f && nlo < 8) off_lo[nlo++] = dlo;
        f = 0;
        for (long j = 0; j < nhi; j++) if (off_hi[j] == dhi) f = 1;
        if (!f && nhi < 8) off_hi[nhi++] = dhi;
        if (dlo == 0 && dhi == 0) exact++;
        if (s1 >= b0[k] && s0 <= b1[k]) contained++;
    }
    printf("   distinct sigma(b1_e) - b0_k : ");
    for (long j = 0; j < nlo; j++) printf("%ld ", off_lo[j]);
    printf("%s\n", nlo >= 8 ? "(>=8 distinct)" : "");
    printf("   distinct sigma(b0_e) - b1_k : ");
    for (long j = 0; j < nhi; j++) printf("%ld ", off_hi[j]);
    printf("%s\n", nhi >= 8 ? "(>=8 distinct)" : "");
    printf("   exact endpoint match: %ld/%ld   reflected interval inside partner: %ld/%ld\n",
           exact, nex, contained, nex);

    if (getenv("R27_DUMP")) {
        printf("#DUMP e b0 b1 len jump off_lo off_hi\n");
        for (long i = 0; i < nex; i++) {
            long k = nex - 1 - i;
            long s0 = ((C - b0[i]) % p + p) % p;
            long s1 = ((C - b1[i]) % p + p) % p;
            printf("#D %ld %ld %ld %ld %d %ld %ld\n", i, b0[i], b1[i],
                   b1[i] - b0[i], jm[i], s1 - b0[k], s0 - b1[k]);
        }
    }

    printf("\n   first 6 excursions: b0, b1, sigma(b1), and the partner's [b0,b1]\n");
    for (long i = 0; i < 6 && i < nex; i++) {
        long k = nex - 1 - i;
        long s1 = ((C - b1[i]) % p + p) % p;
        printf("   e=%-3ld b0=%-9ld b1=%-9ld sigma(b1)=%-9ld  partner e=%-3ld [%ld, %ld]  delta=%ld\n",
               i, b0[i], b1[i], s1, k, b0[k], b1[k], s1 - b0[k]);
    }
    return 0;
}
