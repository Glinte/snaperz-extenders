/* Round-27c: are the 83 dirty excursions of the 56|13 cut exactly the 83 rare
 * large gaps of the exceptional-ECO edge set?
 *
 * Motivation.  r27_reflect showed the clean 56|13 checkpoint section is totally
 * disjoint from its reflection (0/268,128), and that the excursion intervals
 * are reflected only up to an odd offset <= 29, with 76 of 83 partner lengths
 * differing.  So the 84-excursion palindrome is not an interval reflection.
 *
 * But r27b_eco showed the excursions contain almost no exceptional edges, and
 * the multiplicities coincide exactly:
 *
 *     excursion jumps  : 30^2  39^40  14^40  123^1     (83)
 *     Xi weights       : 77^2  379^40 2597^40 254167^1 (83 rare, of 29,265)
 *
 * If excursion e sits in the e-th rare gap, order-preservingly, then the
 * excursion palindrome is a corollary of the Xi-weight palindrome, which the
 * round-26 reversor already proves.  That would close the last L=69 conjecture.
 *
 * usage: r27c_bridge [L] [cut]      (defaults 69 56)
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
    int m = n - 1; while (m >= 0 && !x[m]) m--; return m;
}
static int hgt(const uint8_t *x, int n) { return n - last_pos(x, n); }

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
static int is_exceptional(const uint8_t *x, uint8_t *t, uint8_t *at) {
    int r; unlam(x, L, t, &r);
    int ht = hgt(t, L - 1);
    memcpy(at, t, L - 1); sweep_n(at, L - 1);
    return (hgt(at, L - 1) == ht - 1) && (r == ht);
}
static void build_cycle(void) {
    uint8_t *x = calloc(S, 1), *x0 = calloc(S, 1);
    for (int i = 0; i < S; i++) x0[i] = 1;
    memcpy(x, x0, S);
    cyc = malloc((size_t)(1 << 16) * S); pS = 0;
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
                if (ph >= 0) { is_clean[s] = 1; nclean++;
                    dat[s] = (int32_t)((((long)ph - s) % pS + pS) % pS); }
            }
        }
        sweep_n(x, L);
    }

    /* exceptional-edge positions and their cyclic gaps (the Xi weights) */
    long *g = malloc(nexc * sizeof(long)); long ng = 0;
    for (long s = 0; s < p; s++) if (exc[s]) g[ng++] = s;
    long *w = malloc(ng * sizeof(long));
    for (long i = 0; i < ng; i++) w[i] = (i + 1 < ng ? g[i + 1] : g[0] + p) - g[i];

    /* rare = anything above the four common weights */
    long nrare = 0;
    for (long i = 0; i < ng; i++) if (w[i] > 31) nrare++;
    printf("L=%d cut=%d|%d  p_L=%ld  exceptional edges=%ld  gaps>31: %ld\n",
           L, CUT, S, p, ng, nrare);

    /* excursions */
    long *ct = malloc(nclean * sizeof(long)), m = 0;
    for (long s = 0; s < p; s++) if (is_clean[s]) ct[m++] = s;
    long nex = 0;
    for (long i = 1; i < m; i++) if (dat[ct[i]] != dat[ct[i - 1]]) nex++;
    long *b0 = malloc(nex * sizeof(long)), *b1 = malloc(nex * sizeof(long));
    int *jm = malloc(nex * sizeof(int));
    long e = 0;
    for (long i = 1; i < m; i++)
        if (dat[ct[i]] != dat[ct[i - 1]]) {
            b0[e] = ct[i - 1]; b1[e] = ct[i];
            jm[e] = (int)(((dat[ct[i]] - dat[ct[i - 1]]) % pS + pS) % pS);
            e++;
        }
    printf("excursions=%ld   rare gaps=%ld   %s\n", nex, nrare,
           nex == nrare ? "COUNTS MATCH" : "COUNTS DIFFER");

    /* Which gap contains each excursion?  r27b_eco measured b1 to be within 71
     * of an exceptional edge always, while b0 can be 254,151 away, so b1 is the
     * anchor with meaning; ANCHOR=b0 re-runs the weaker version. */
    const char *anch = getenv("ANCHOR");
    int use_b0 = anch && !strcmp(anch, "b0");
    long *gidx = malloc(nex * sizeof(long));
    for (long i = 0; i < nex; i++) {
        long a = use_b0 ? b0[i] : b1[i];
        long lo = 0, hi = ng - 1, k = ng - 1;
        while (lo <= hi) { long mid = (lo + hi) / 2;
            if (g[mid] <= a) { k = mid; lo = mid + 1; } else hi = mid - 1; }
        gidx[i] = k;
    }
    printf("anchor: %s\n", use_b0 ? "b0" : "b1");
    long rare_hit = 0, distinct = 1, monotone = 1;
    for (long i = 0; i < nex; i++) {
        if (w[gidx[i]] > 31) rare_hit++;
        if (i && gidx[i] <= gidx[i - 1]) monotone = 0;
        if (i && gidx[i] == gidx[i - 1]) distinct = 0;
    }
    printf("excursion b0 lands in a rare gap: %ld/%ld   distinct gaps: %s   order-preserving: %s\n",
           rare_hit, nex, distinct ? "yes" : "NO", monotone ? "yes" : "NO");

    /* the enclosing weight, and the weight -> jump map */
    printf("\n  enclosing gap weight -> jump (single-valued?):\n");
    long uw[64]; int uj[64]; long uc[64], nu = 0; int multi = 0;
    for (long i = 0; i < nex; i++) {
        long ww = w[gidx[i]];
        int f = -1;
        for (long j = 0; j < nu; j++) if (uw[j] == ww) f = (int)j;
        if (f < 0 && nu < 64) { f = (int)nu++; uw[f] = ww; uj[f] = jm[i]; uc[f] = 0; }
        if (f >= 0) { uc[f]++; if (uj[f] != jm[i]) multi = 1; }
    }
    for (long j = 0; j < nu; j++)
        printf("    weight=%-8ld count=%-4ld jump=%d\n", uw[j], uc[j], uj[j]);
    printf("    single-valued: %s\n", multi ? "NO" : "YES");

    /* Do not assume a rooting.  The exceptional set is sigma-invariant, so
     * gap i = [g_i, g_{i+1}) reflects onto the gap starting at C - g_{i+1}.
     * Build that involution pi directly and check w_{pi(i)} = w_i. */
    long C = ((p - 2 * (long)L + 2) % p + p) % p;
    long *pos = malloc(p * sizeof(long));
    for (long i = 0; i < p; i++) pos[i] = -1;
    for (long i = 0; i < ng; i++) pos[g[i]] = i;

    long *pi = malloc(ng * sizeof(long));
    long pi_ok = 0, wpal_true = 0;
    for (long i = 0; i < ng; i++) {
        long nxt = g[(i + 1) % ng];
        long img = ((C - nxt) % p + p) % p;
        pi[i] = pos[img];
        if (pi[i] >= 0) { pi_ok++; if (w[pi[i]] == w[i]) wpal_true++; }
    }
    printf("\n  reflection on gaps: pi well-defined %ld/%ld, w_{pi(i)} = w_i on %ld/%ld\n",
           pi_ok, ng, wpal_true, pi_ok);

    /* induced involution on the 83 rare gaps, in cyclic order */
    long *rare = malloc(ng * sizeof(long)), nr = 0;
    long *rank = malloc(ng * sizeof(long));
    for (long i = 0; i < ng; i++) rank[i] = -1;
    for (long i = 0; i < ng; i++) if (w[i] > 31) { rank[i] = nr; rare[nr++] = i; }

    long affine = 0, cst = -1; int consistent = 1;
    for (long r = 0; r < nr; r++) {
        long j = pi[rare[r]];
        if (j < 0 || rank[j] < 0) continue;
        long c = (r + rank[j]) % nr;
        if (cst < 0) cst = c;
        else if (c != cst) consistent = 0;
        affine++;
    }
    printf("  rare gaps: %ld, reflection sends rank r -> %s (const %ld), on %ld of them\n",
           nr, consistent ? "const - r" : "NO affine form", cst, affine);

    /* does the excursion pairing e <-> nex-1-e agree with the rare-gap one? */
    long agree = 0;
    for (long i = 0; i < nex; i++) {
        long r = rank[gidx[i]];
        if (r < 0) continue;
        long rp = rank[gidx[nex - 1 - i]];
        if (rp < 0) continue;
        if ((r + rp) % nr == cst) agree++;
    }
    printf("  excursion pairing agrees with the rare-gap reflection: %ld/%ld\n", agree, nex);

    long jpal = 0;
    for (long i = 0; i < nex; i++) if (jm[i] == jm[nex - 1 - i]) jpal++;
    printf("  jump word palindromic: %ld/%ld\n", jpal, nex);

    /* Better anchor: the rare gap the excursion interval actually overlaps.
     * If that is unique per excursion the correspondence is a genuine
     * bijection, and the excursion palindrome descends from the gap one. */
    long uniq = 0, none = 0, many = 0;
    long *ov = malloc(nex * sizeof(long));
    for (long i = 0; i < nex; i++) {
        long hits = 0, last = -1;
        for (long r = 0; r < nr; r++) {
            long j = rare[r], s0 = g[j], s1 = s0 + w[j];
            if (s0 < b1[i] && s1 > b0[i]) { hits++; last = r; }
        }
        ov[i] = hits == 1 ? last : -1;
        if (hits == 1) uniq++; else if (!hits) none++; else many++;
    }
    printf("\n  rare gaps overlapping each excursion: unique %ld, none %ld, several %ld (of %ld)\n",
           uniq, none, many, nex);
    if (uniq == nex) {
        long mono = 1, ag = 0;
        for (long i = 1; i < nex; i++) if (ov[i] <= ov[i - 1]) mono = 0;
        for (long i = 0; i < nex; i++)
            if ((ov[i] + ov[nex - 1 - i]) % nr == cst) ag++;
        printf("  order-preserving: %s   pairing agrees with gap reflection: %ld/%ld\n",
               mono ? "yes" : "NO", ag, nex);
    }
    return 0;
}
