/* Round-14 audit, part F: the 56+s tail phase-lock table and the injected root q. */
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

static int cyc[2000][40], nc;
static void build(int s) {
    int y[40], y0[40];
    for (int i = 0; i < s; i++) y[i] = y0[i] = 1;
    nc = 0; memcpy(cyc[nc++], y, s * sizeof(int));
    for (;;) { Astep(y, s); if (!memcmp(y, y0, s * sizeof(int))) break; memcpy(cyc[nc++], y, s * sizeof(int)); }
}
static int phase(const int *t, int s) {
    for (int i = 0; i < nc; i++) if (!memcmp(cyc[i], t, s * sizeof(int))) return i;
    return -1;
}

int main(void) {
    for (int s = 9; s <= 17; s += 2) {
        int L = 56 + s;
        build(s);
        int ps = nc;
        int x[MAXL], x0[MAXL];
        for (int i = 0; i < L; i++) x[i] = x0[i] = 1;
        long p = 0, k = 0;
        int nclean = 0, noncyc = 0, lockfail = 0, c_s = -999, first = 1;
        int phi780 = -1, q = -1, T780 = -1;
        for (;;) {
            Astep(x, L); p++;
            if (x[0] == 1) {
                k++;
                if (k % 112 == 0) {
                    long j = k / 112;
                    if (j > 780) break;
                    int sum = 0; for (int i = 0; i < 56; i++) sum += x[i];
                    if (sum == 56) {
                        nclean++;
                        int ph = phase(x + 56, s);
                        if (ph < 0) noncyc++;
                        else {
                            int c = ((ph - (int)(p % ps)) % ps + ps) % ps;
                            if (first) { c_s = c; first = 0; }
                            else if (c != c_s) lockfail++;
                            if (j == 780) { phi780 = ph; T780 = (int)p;
                                            int t[40]; memcpy(t, cyc[ph], s * sizeof(int));
                                            for (int u = 0; u < 293; u++) Astep(t, s);
                                            q = t[0]; }
                        }
                    }
                }
            }
            if (!memcmp(x, x0, L * sizeof(int))) break;
        }
        printf("s=%2d  L=%2d  p_s=%3d  clean cuts (j<=780) %3d  tails off cycle %d  "
               "phase-lock phi=T+c_s fails %d  c_s=%d  phi_780=%d  T_780=%d  "
               "q=(A^293 y)_0=%d (%s)\n",
               s, L, ps, nclean, noncyc, lockfail, c_s, phi780, T780, q, q % 2 ? "odd" : "even");
    }
    return 0;
}
