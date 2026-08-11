/* Round-16 audit: the boundary-event cocycle of the clean 56|13 cut at L=69.
 *
 * Rebuilds the canonical A_69 orbit from the gate rule, records the boundary
 * symbol seen by cell 55 on every sweep, and reconstructs the claimed
 * constant-d runs, dirty excursions, event-word alphabet, palindromic jump
 * word, holonomy, and phase-reflection law.
 *
 * usage: r16_cocycle [L] [cut]        (defaults 69 56)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int L, CUT, S;                 /* prefix length CUT, tail length S */

static uint8_t *cyc;                  /* canonical A_S cycle, pS states */
static int pS;

static void sweep_n(uint8_t *x, int n) {
    for (int i = 0; i < n - 1; i++) {
        uint8_t a = x[i];
        if (!a) continue;
        if (a == 1) { x[i] = 1 + x[i + 1]; x[i + 1] = 0; }
        else { x[i] = a - 1; x[i + 1] += 1; }
    }
}

/* one L-sweep, returning the boundary symbol seen at gate CUT-1 (0/1/2) */
static int sweep_boundary(uint8_t *x) {
    int sym = -1;
    for (int i = 0; i < L - 1; i++) {
        if (i == CUT - 1) sym = x[i] >= 2 ? 2 : x[i];
        uint8_t a = x[i];
        if (!a) continue;
        if (a == 1) { x[i] = 1 + x[i + 1]; x[i + 1] = 0; }
        else { x[i] = a - 1; x[i + 1] += 1; }
    }
    return sym;
}

static void build_cycle(void) {
    uint8_t *x = calloc(S, 1), *x0 = calloc(S, 1);
    for (int i = 0; i < S; i++) x0[i] = 1;
    memcpy(x, x0, S);
    int cap = 1 << 16;
    cyc = malloc((size_t)cap * S);
    pS = 0;
    do {
        memcpy(cyc + (size_t)pS * S, x, S);
        pS++;
        sweep_n(x, S);
    } while (memcmp(x, x0, S));
    free(x); free(x0);
}

static int phase_of(const uint8_t *tail) {
    for (int k = 0; k < pS; k++)
        if (!memcmp(cyc + (size_t)k * S, tail, S)) return k;
    return -1;
}

static int zdepth(const uint8_t *x, int n) {
    for (int i = 0; i < n; i++) if (!x[i]) return i;
    return n;
}

int main(int argc, char **argv) {
    L = argc > 1 ? atoi(argv[1]) : 69;
    CUT = argc > 2 ? atoi(argv[2]) : 56;
    S = L - CUT;
    build_cycle();
    printf("L=%d cut=%d|%d  p_%d=%d\n", L, CUT, S, S, pS);

    uint8_t *x = calloc(L, 1), *x0 = calloc(L, 1);
    for (int i = 0; i < L; i++) x0[i] = 1;
    memcpy(x, x0, L);

    long cap = 1 << 22, p = 0;
    uint8_t *sym = malloc(cap);
    long *cp_t = malloc(cap * sizeof(long));     /* checkpoint times */
    int32_t *cp_ph = malloc(cap * sizeof(int32_t)); /* tail phase, -1 = none */
    int32_t *cp_z = malloc(cap * sizeof(int32_t));
    uint8_t *cp_clean = malloc(cap);
    long nc = 0;

    /* index 0 is E_L itself: a checkpoint at t = 0 */
    for (;;) {
        if (x[0] == 1) {
            int s = 0, clean = 1;
            for (int i = 0; i < CUT; i++) { s += x[i]; if (s < i + 1) clean = 0; }
            if (s != CUT) clean = 0;
            cp_clean[nc] = (uint8_t)clean;
            cp_t[nc] = p;
            cp_z[nc] = zdepth(x, L);
            cp_ph[nc] = clean ? phase_of(x + CUT) : -1;
            nc++;
        }
        sym[p] = (uint8_t)sweep_boundary(x);
        p++;
        if (!memcmp(x, x0, L)) break;
    }
    printf("p_%d=%ld checkpoints=%ld\n", L, p, nc);

    long clean_n = 0, oncyc = 0;
    for (long k = 0; k < nc; k++) {
        if (cp_clean[k]) clean_n++;
        if (cp_ph[k] >= 0) oncyc++;
    }
    printf("clean %d-cuts=%ld  tail on canonical A_%d cycle=%ld  off-cycle=%ld\n",
           CUT, clean_n, S, oncyc, clean_n - oncyc);

    /* constant-d runs over the clean on-cycle checkpoints, in time order */
    long *idx = malloc(nc * sizeof(long));
    long m = 0;
    for (long k = 0; k < nc; k++) if (cp_ph[k] >= 0) idx[m++] = k;

    int *dval = malloc(m * sizeof(int));
    for (long i = 0; i < m; i++) {
        long k = idx[i];
        dval[i] = (int)(((cp_ph[k] - cp_t[k]) % pS + pS) % pS);
    }
    long runs = 1;
    for (long i = 1; i < m; i++) if (dval[i] != dval[i - 1]) runs++;
    printf("constant-d runs=%ld  excursions=%ld  d_first=%d d_last=%d\n",
           runs, runs - 1, dval[0], dval[m - 1]);

    /* event words of each excursion: the non-Z boundary symbols strictly
       between the last checkpoint of one run and the first of the next */
    long nex = runs - 1;
    long *elen = malloc(nex * sizeof(long));
    uint64_t *ehash = malloc(nex * sizeof(uint64_t));
    int *ejump = malloc(nex * sizeof(int));
    long *ends = malloc(nex * sizeof(long));
    long e = 0;
    for (long i = 1; i < m; i++) {
        if (dval[i] == dval[i - 1]) continue;
        long t0 = cp_t[idx[i - 1]], t1 = cp_t[idx[i]];
        uint64_t hsh = 1469598103934665603ULL;
        long len = 0;
        for (long t = t0; t < t1; t++)
            if (sym[t]) { hsh = (hsh ^ sym[t]) * 1099511628211ULL; len++; }
        elen[e] = len; ehash[e] = hsh;
        ejump[e] = ((dval[i] - dval[i - 1]) % pS + pS) % pS;
        ends[e] = i;
        e++;
    }

    /* distinct event words */
    printf("\ndistinct event words (length, count, jump):\n");
    uint64_t seen[64]; long slen[64], scnt[64]; int sj[64], ns = 0;
    for (long i = 0; i < nex; i++) {
        int f = -1;
        for (int j = 0; j < ns; j++) if (seen[j] == ehash[i] && slen[j] == elen[i]) f = j;
        if (f < 0 && ns < 64) { f = ns++; seen[f] = ehash[i]; slen[f] = elen[i]; scnt[f] = 0; sj[f] = ejump[i]; }
        if (f >= 0) { scnt[f]++; if (sj[f] != ejump[i]) sj[f] = -1; }
    }
    for (int j = 0; j < ns; j++)
        printf("  len=%-8ld count=%-5ld jump=%d\n", slen[j], scnt[j], sj[j]);
    if (ns >= 64) printf("  (more than 64 distinct words)\n");

    /* jump word and its palindromy */
    int palj = 1;
    long S_tot = 0;
    for (long i = 0; i < nex; i++) {
        S_tot += ejump[i];
        if (ejump[i] != ejump[nex - 1 - i]) palj = 0;
    }
    printf("\njump word palindromic=%d  sum=%ld  sum mod %d=%ld\n",
           palj, S_tot, pS, S_tot % pS);
    printf("jump word: ");
    for (long i = 0; i < nex && i < 100; i++) printf("%d ", ejump[i]);
    printf("%s\n", nex > 100 ? "..." : "");

    /* defect reflection d_i + d_{m-1-i} and the phase-reflection law */
    long refl_bad = 0, refl_n = 0;
    int hol = (dval[0] + dval[m - 1]) % pS;
    for (long i = 0; i < m; i++) {
        long j = m - 1 - i;
        refl_n++;
        if ((dval[i] + dval[j]) % pS != hol) refl_bad++;
    }
    printf("\nrun-index defect reflection: pairs=%ld failures=%ld  holonomy=%d\n",
           refl_n, refl_bad, hol);

    /* phi_k + phi_{N-k} = z(P_k) + C mod pS over clean checkpoint pairs */
    long ok = 0, bad = 0; int Cconst = -1;
    for (long k = 1; k < nc; k++) {
        long k2 = nc - k;
        if (cp_ph[k] < 0 || cp_ph[k2] < 0) continue;
        int c = ((cp_ph[k] + cp_ph[k2] - cp_z[k]) % pS + pS) % pS;
        if (Cconst < 0) Cconst = c;
        if (c == Cconst) ok++; else bad++;
    }
    printf("phase reflection phi_k+phi_{N-k} = z(x_k)+C: C=%d ok=%ld bad=%ld\n",
           Cconst, ok, bad);
    printf("check: (p_L - L) mod %d = %ld,  + holonomy %d = %d\n",
           pS, ((p - L) % pS + pS) % pS, hol, (int)((((p - L) % pS + pS) % pS + hol) % pS));
    return 0;
}
