/* Round-29: is the excursion <-> rare-Xi-gap correspondence an L=69 accident,
 * or the general mechanism at 56|s cuts?
 *
 * r27c_bridge established at L=69, 56|13: the 83 dirty excursions sit inside
 * the 83 rare exceptional-edge gaps, order-preservingly, with the enclosing
 * weight determining the jump.  Round 24's cut audits give defect-run counts
 * at other lengths (L=64 56|8: 937 runs; L=68 56|12: 48; L=70 56|14: 230), so
 * the same question is decidable there: does #excursions equal #gaps above a
 * threshold, and does each excursion overlap exactly one rare gap?
 *
 * Memory-lean rewrite of r27c: stores only event lists (exceptional-edge
 * times, clean-checkpoint times + defects), never a per-sweep array, so a
 * 73M-sweep orbit fits in a 3GB machine.
 *
 * usage: r29_cuts L cut [rare_threshold]      (threshold default: report only)
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
    cyc = malloc((size_t)(1 << 16) * S);
    pS = 0;
    do { memcpy(cyc + (size_t)pS * S, x, S); pS++; sweep_n(x, S); }
    while (memcmp(x, x0, S));
    free(x); free(x0);
}

/* tiny open-addressed hash: tail bytes -> phase */
static uint32_t *htab; static int hcap;
static uint64_t thash(const uint8_t *t) {
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < S; i++) h = (h ^ t[i]) * 1099511628211ULL;
    return h;
}
static void hbuild(void) {
    hcap = 4; while (hcap < 4 * pS) hcap <<= 1;
    htab = malloc((size_t)hcap * 4);
    memset(htab, 0xff, (size_t)hcap * 4);
    for (int k = 0; k < pS; k++) {
        uint64_t h = thash(cyc + (size_t)k * S) & (hcap - 1);
        while (htab[h] != 0xffffffffu) h = (h + 1) & (hcap - 1);
        htab[h] = (uint32_t)k;
    }
}
static int phase_of(const uint8_t *t) {
    uint64_t h = thash(t) & (hcap - 1);
    while (htab[h] != 0xffffffffu) {
        uint32_t k = htab[h];
        if (!memcmp(cyc + (size_t)k * S, t, S)) return (int)k;
        h = (h + 1) & (hcap - 1);
    }
    return -1;
}

static long bfind(const long *a, long n, long v) {   /* index of v, or -1 */
    long lo = 0, hi = n - 1;
    while (lo <= hi) {
        long mid = (lo + hi) / 2;
        if (a[mid] == v) return mid;
        if (a[mid] < v) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

int main(int argc, char **argv) {
    L = argc > 1 ? atoi(argv[1]) : 69;
    CUT = argc > 2 ? atoi(argv[2]) : 56;
    long rareT = argc > 3 ? atol(argv[3]) : 0;
    S = L - CUT;
    build_cycle(); hbuild();

    uint8_t *x = calloc(L, 1), *x0 = calloc(L, 1);
    uint8_t *t = calloc(L, 1), *at = calloc(L, 1);
    for (int i = 0; i < L; i++) x0[i] = 1;
    memcpy(x, x0, L);

    long gcap = 1 << 20, ccap = 1 << 20;
    long *g = malloc(gcap * 8);
    long *ct = malloc(ccap * 8);
    int32_t *cd = malloc(ccap * 4);
    int64_t *ck = malloc(ccap * 8);   /* global checkpoint index of each clean one */
    long ng = 0, m = 0, ncp = 0, nclean_off = 0;

    long p = 0;
    do {
        if (is_exceptional(x, t, at)) {
            if (ng == gcap) { gcap <<= 1; g = realloc(g, gcap * 8); }
            g[ng++] = p;
        }
        if (x[0] == 1) {
            ncp++;
            int sm = 0, clean = 1;
            for (int i = 0; i < CUT; i++) { sm += x[i]; if (sm < i + 1) clean = 0; }
            if (sm != CUT) clean = 0;
            if (clean) {
                int ph = phase_of(x + CUT);
                if (ph >= 0) {
                    if (m == ccap) { ccap <<= 1; ct = realloc(ct, ccap * 8); cd = realloc(cd, ccap * 4); ck = realloc(ck, ccap * 8); }
                    ct[m] = p;
                    cd[m] = (int32_t)((((long)ph - p) % pS + pS) % pS);
                    ck[m] = ncp - 1;
                    m++;
                } else nclean_off++;
            }
        }
        sweep_n(x, L);
        p++;
    } while (memcmp(x, x0, L));

    printf("L=%d cut=%d|%d  p_L=%ld  p_%d=%d\n", L, CUT, S, p, S, pS);
    printf("checkpoints=%ld  clean on-cycle=%ld (off-cycle clean=%ld)  exceptional edges=%ld\n",
           ncp, m, nclean_off, ng);

    /* sanity: the edge set must be sigma-invariant (reversor theorem) */
    long C = ((p - 2 * (long)L + 2) % p + p) % p;
    long sym = 0;
    for (long i = 0; i < ng; i++)
        if (bfind(g, ng, ((C - g[i]) % p + p) % p) >= 0) sym++;
    printf("edge set sigma-invariant: %ld/%ld  (centre %ld)\n", sym, ng, C);

    /* gap weights and histogram */
    long *w = malloc(ng * 8);
    for (long i = 0; i < ng; i++) w[i] = (i + 1 < ng ? g[i + 1] : g[0] + p) - g[i];
    long uw[512], uc[512], nu = 0;
    for (long i = 0; i < ng; i++) {
        long f = -1;
        for (long j = 0; j < nu; j++) if (uw[j] == w[i]) f = j;
        if (f < 0 && nu < 512) { f = nu++; uw[f] = w[i]; uc[f] = 0; }
        if (f >= 0) uc[f]++;
    }
    for (long a = 0; a < nu; a++)          /* sort by weight */
        for (long b = a + 1; b < nu; b++)
            if (uw[b] < uw[a]) {
                long tt = uw[a]; uw[a] = uw[b]; uw[b] = tt;
                tt = uc[a]; uc[a] = uc[b]; uc[b] = tt;
            }
    printf("gap weights: ");
    for (long j = 0; j < nu; j++) printf("%ld^%ld ", uw[j], uc[j]);
    printf("%s\n", nu >= 512 ? "(truncated)" : "");

    /* excursions */
    long nex = 0;
    for (long i = 1; i < m; i++) if (cd[i] != cd[i - 1]) nex++;
    printf("excursions=%ld\n", nex);
    printf("gaps above w: ");
    long cum = 0;
    for (long j = nu - 1; j >= 0 && nu - j <= 12; j--) {
        cum += uc[j];
        printf("[>%ld]:%ld ", j ? uw[j - 1] : 0, cum);
    }
    printf("\n");
    if (!rareT || !nex) return 0;

    long *b0 = malloc(nex * 8), *b1 = malloc(nex * 8);
    int *jm = malloc(nex * 4);
    long e = 0;
    for (long i = 1; i < m; i++)
        if (cd[i] != cd[i - 1]) {
            b0[e] = ct[i - 1]; b1[e] = ct[i];
            jm[e] = (int)(((cd[i] - cd[i - 1]) % pS + pS) % pS);
            e++;
        }

    /* rare gaps and the overlap correspondence */
    long nr = 0;
    for (long i = 0; i < ng; i++) if (w[i] > rareT) nr++;
    printf("\nthreshold %ld: rare gaps=%ld vs excursions=%ld  %s\n",
           rareT, nr, nex, nr == nex ? "COUNTS MATCH" : "COUNTS DIFFER");

    long uniq = 0, none = 0, many = 0, mono = 1, prev = -1;
    int wj_multi = 0;
    long wjw[64]; int wjj[64]; long wjc[64], nwj = 0;
    for (long i = 0; i < nex; i++) {
        long hits = 0, hitg = -1;
        /* gaps intersecting [b0, b1): binary search the first g >= b0, then
           also the gap containing b0 */
        long lo = 0, hi = ng - 1, k0 = ng - 1;
        while (lo <= hi) { long mid = (lo + hi) / 2;
            if (g[mid] <= b0[i]) { k0 = mid; lo = mid + 1; } else hi = mid - 1; }
        for (long k = k0; k < ng && g[k] < b1[i]; k++)
            if (w[k] > rareT && g[k] + w[k] > b0[i]) { hits++; hitg = k; }
        if (hits == 1) {
            uniq++;
            if (hitg <= prev) mono = 0;
            prev = hitg;
            long ww = w[hitg];
            long f = -1;
            for (long j = 0; j < nwj; j++) if (wjw[j] == ww) f = j;
            if (f < 0 && nwj < 64) { f = nwj++; wjw[f] = ww; wjj[f] = jm[i]; wjc[f] = 0; }
            if (f >= 0) { wjc[f]++; if (wjj[f] != jm[i]) wj_multi = 1; }
        } else if (!hits) none++;
        else many++;
    }
    printf("excursion overlaps exactly one rare gap: %ld/%ld (none %ld, several %ld)  order-preserving: %s\n",
           uniq, nex, none, many, mono ? "yes" : "NO");
    printf("weight -> jump: ");
    for (long j = 0; j < nwj; j++) printf("%ld->%d(x%ld) ", wjw[j], wjj[j], wjc[j]);
    printf("  single-valued: %s\n", wj_multi ? "NO" : "YES");

    long jpal = 0;
    for (long i = 0; i < nex; i++) if (jm[i] == jm[nex - 1 - i]) jpal++;
    printf("jump word palindromic: %ld/%ld\n", jpal, nex);

    /* The proposed two-line derivation: reflection law + time identity give
     *     d_k + d_{N-k} = C - (p_L - L)   (mod pS),  a single constant,
     * with k the GLOBAL checkpoint index, N = ncp.  Measure it directly. */
    long hist[512]; memset(hist, 0, sizeof hist);
    long paired = 0, orphan = 0;
    for (long i = 0; i < m; i++) {
        long want = (ncp - ck[i]) % ncp;
        /* binary search ck for want */
        long lo = 0, hi = m - 1, j = -1;
        while (lo <= hi) { long mid = (lo + hi) / 2;
            if (ck[mid] == want) { j = mid; break; }
            if (ck[mid] < want) lo = mid + 1; else hi = mid - 1; }
        if (j < 0) { orphan++; continue; }
        paired++;
        hist[(cd[i] + cd[j]) % pS]++;
    }
    long best = 0, bestk = -1, nkv = 0;
    for (int k = 0; k < pS && k < 512; k++)
        if (hist[k]) { nkv++; if (hist[k] > best) { best = hist[k]; bestk = k; } }
    printf("\ndefect-sum law d_k + d_(N-k) mod %d over %ld clean pairs (%ld orphans):\n",
           pS, paired, orphan);
    printf("  distinct values: %ld   mode %ld with %ld (%.3f%%)   predicted C-(p_L-L) mod %d: ",
           nkv, bestk, best, 100.0 * best / (paired ? paired : 1), pS);
    printf("(fill in C by hand)\n");
    if (nkv <= 8) {
        printf("  values: ");
        for (int k = 0; k < pS && k < 512; k++)
            if (hist[k]) printf("%d:%ld ", k, hist[k]);
        printf("\n");
    }
    return 0;
}
