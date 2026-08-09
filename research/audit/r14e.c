/* Round-14 audit, part E: the controlled 129-phase experiment at L=69,
   and the common 56-cell prefix across odd neighbours at outer step 780. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXL 128

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

/* run to outer step J (J*112 checkpoints) from E_L, copying the state out */
static int section_at(int L, long J, int *out) {
    int x[MAXL], x0[MAXL];
    for (int i = 0; i < L; i++) x[i] = x0[i] = 1;
    long k = 0;
    if (J == 0) { memcpy(out, x, L * sizeof(int)); return 1; }
    for (;;) {
        Astep(x, L);
        if (x[0] == 1) {
            k++;
            if (k == J * 112) { memcpy(out, x, L * sizeof(int)); return 1; }
        }
        if (!memcmp(x, x0, L * sizeof(int))) return 0;   /* closed early */
    }
}

int main(void) {
    /* ---- 1. common 56-cell prefix at outer step 780 across odd L ---- */
    printf("--- outer step 780, odd neighbours ---\n");
    static int ref[MAXL];
    int haveref = 0;
    for (int L = 65; L <= 83; L += 2) {
        static int st[MAXL];
        if (!section_at(L, 780, st)) { printf("L=%d closed before outer step 780\n", L); continue; }
        if (!haveref) { memcpy(ref, st, 56 * sizeof(int)); haveref = 1; }
        int same = !memcmp(ref, st, 56 * sizeof(int));
        int tailsum = 0; for (int i = 56; i < L; i++) tailsum += st[i];
        printf("L=%2d  first 56 coords match L=65: %s   tail root q = x_56 = %d (%s)   ell = %d\n",
               L, same ? "yes" : "NO", st[56], st[56] % 2 ? "odd" : "even", ell(st, L));
        (void)tailsum;
    }

    /* ---- 2. the canonical A_13 cycle ---- */
    static int cyc[200][13];
    int nc = 0;
    { int y[13], y0[13];
      for (int i = 0; i < 13; i++) y[i] = y0[i] = 1;
      memcpy(cyc[nc++], y, sizeof y);
      for (;;) { Astep(y, 13); if (!memcmp(y, y0, sizeof y)) break; memcpy(cyc[nc++], y, sizeof y); } }
    printf("\nA_13 canonical cycle length %d\n", nc);

    /* ---- 3. 129-phase experiment against the real P_56 at outer step 780 ---- */
    static int X780[MAXL];
    if (!section_at(69, 780, X780)) { printf("L=69 closed before 780\n"); return 1; }
    printf("real L=69 tail at j=780: ");
    for (int i = 56; i < 69; i++) printf("%d ", X780[i]);
    int realphase = -1;
    for (int p = 0; p < nc; p++) if (!memcmp(cyc[p], X780 + 56, 13 * sizeof(int))) realphase = p;
    printf("  -> phase %d\n\n", realphase);

    printf("--- 129 controlled trajectories (P_56 fixed, tail = every A_13 phase) ---\n");
    int surv834 = 0, surv1614 = 0, reachE = 0;
    int surv834_list[200], nsl = 0, reach_list[200], nrl = 0;
    for (int p = 0; p < nc; p++) {
        int x[MAXL], E69[MAXL];
        for (int i = 0; i < 69; i++) E69[i] = 1;
        memcpy(x, X780, 56 * sizeof(int));
        memcpy(x + 56, cyc[p], 13 * sizeof(int));
        int ok834 = 1, ok1614 = 1, hitE = 0; long hitJ = -1;
        long cps = 0, sweeps = 0;
        for (long j = 1; j <= 1614; j++) {
            /* one outer return = 112 checkpoints */
            for (int c = 0; c < 112; c++) {
                for (;;) {
                    Astep(x, 69); sweeps++;
                    if (!memcmp(x, E69, 69 * sizeof(int)) && !hitE) { hitE = 1; hitJ = j; }
                    if (x[0] == 1) { cps++; break; }
                    if (sweeps > 20000000L) break;
                }
            }
            int l = ell(x, 69);
            if (l < 56) { if (j <= 834) ok834 = 0; ok1614 = 0; }
        }
        if (ok834) { surv834++; surv834_list[nsl++] = p; }
        if (ok1614) surv1614++;
        if (hitE) { reachE++; reach_list[nrl++] = p; if (nrl < 8) printf("  phase %3d reaches E_69 at outer return %ld\n", p, hitJ); }
    }
    printf("phases surviving ell>=56 for the first 834 outer returns: %d  ->", surv834);
    for (int i = 0; i < nsl; i++) printf(" %d", surv834_list[i]);
    printf("   (claim {19,21,23,25})\n");
    printf("phases surviving all 1614 outer returns: %d   (claim 3: {19,21,23})\n", surv1614);
    printf("phases reaching E_69 within 1614 outer returns: %d ->", reachE);
    for (int i = 0; i < nrl; i++) printf(" %d", reach_list[i]);
    printf("   (claim exactly {23})\n");
    return 0;
}
