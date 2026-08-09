/* Round-14 audit, part A:
     - L=69 reflection identity  t_k + t_{N-k} = p_L - L + z(x_k)
     - the 56+13 clean-section census and the phase-reflection law
     - the universal terminal cone and its "L=18 gap word" increments
   Built only from the gate rule (same reconstruction as core.py).            */
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

static int zdepth(const int *x, int n) {
    for (int k = 0; k < n; k++) if (x[k] == 0) return k;
    return n;
}

/* ---------- canonical A_18 checkpoint gap word ---------- */
static int g18[256], ng18;
static void build_g18(void) {
    int L = 18, x[18], x0[18];
    for (int i = 0; i < L; i++) x[i] = x0[i] = 1;
    long p = 0, last = 0;
    ng18 = 0;
    for (;;) {
        Astep(x, L); p++;
        if (x[0] == 1) { g18[ng18++] = (int)(p - last); last = p; }
        if (!memcmp(x, x0, sizeof x)) break;
    }
    printf("[g18] period %ld, %d gaps, sum %d\n", p, ng18, (int)p);
}

/* ---------- canonical A_13 cycle, for tail phases ---------- */
static int cyc13[200][13], ncyc13;
static void build_cyc13(void) {
    int x[13], x0[13];
    for (int i = 0; i < 13; i++) x[i] = x0[i] = 1;
    ncyc13 = 0;
    memcpy(cyc13[ncyc13++], x, sizeof x);
    for (;;) {
        Astep(x, 13);
        if (!memcmp(x, x0, sizeof x)) break;
        memcpy(cyc13[ncyc13++], x, sizeof x);
    }
    printf("[A13] canonical cycle length %d\n", ncyc13);
}
static int phase13(const int *t) {
    for (int i = 0; i < ncyc13; i++) if (!memcmp(cyc13[i], t, 13 * sizeof(int))) return i;
    return -1;
}

/* ---------- L = 69 ---------- */
#define NCP 268200
static long tk[NCP + 2];
static int  zk[NCP + 2];
static int  sect[2400][69];

static void run69(void) {
    const int L = 69;
    int x[69], x0[69];
    for (int i = 0; i < L; i++) x[i] = x0[i] = 1;

    long p = 0, k = 0;
    tk[0] = 0; zk[0] = L;                 /* x_0 = E_L, z(E_L) = L */
    int nsec = 0;
    memcpy(sect[nsec++], x, sizeof x);
    for (;;) {
        Astep(x, L); p++;
        if (x[0] == 1) {
            k++;
            tk[k] = p;
            zk[k] = zdepth(x, L);
            if (k % 112 == 0 && !memcmp(x, x0, sizeof x) == 0) { /* placeholder */ }
            if (k % 112 == 0 && k / 112 < 2400) memcpy(sect[k / 112], x, sizeof x);
        }
        if (!memcmp(x, x0, sizeof x)) break;
    }
    long N = k;
    printf("[69] p_69=%ld  B_69=%ld  checkpoints N=%ld  Q_69=%ld  sections=%ld\n",
           p, p - 68, N, N - 1, N / 112);

    /* reflection identity  t_k + t_{N-k} = p - L + z(x_k) */
    long bad = 0, firstbad = -1;
    for (long j = 0; j <= N; j++) {
        long lhs = tk[j] + tk[N - j];
        long rhs = p - L + zk[j];
        if (lhs != rhs) { if (firstbad < 0) firstbad = j; bad++; }
    }
    printf("[69] reflection identity t_k + t_{N-k} = p-L+z(x_k): %s (%ld failures over %ld indices)\n",
           bad ? "FAIL" : "OK", bad, N + 1);

    /* clean 56-cuts, tail phases, z of the prefix */
    int nsecs = (int)(N / 112);
    static int clean[2400], phi[2400], zp[2400];
    int nclean = 0, offcycle = 0;
    for (int j = 0; j < nsecs; j++) {
        int s = 0;
        for (int i = 0; i < 56; i++) s += sect[j][i];
        if (s != 56) { clean[j] = 0; continue; }
        int ph = phase13(sect[j] + 56);
        clean[j] = 1; nclean++;
        if (ph < 0) { offcycle++; phi[j] = -1; }
        else phi[j] = ph;
        zp[j] = zdepth(sect[j], 56);
    }
    printf("[69] clean cuts at 56: %d of %d   (claim 1603)  tails off the canonical A13 cycle: %d (claim 0)\n",
           nclean, nsecs, offcycle);

    /* phase reflection law */
    int paired = 0, nontrivial = 0, zfail = 0, lawfail = 0;
    for (int j = 0; j < nsecs; j++) {
        if (!clean[j]) continue;
        int m = (nsecs - j) % nsecs;
        if (!clean[m]) continue;
        paired++;
        if (j == 0) continue;
        nontrivial++;
        if (zp[m] != zp[j]) zfail++;
        if (((phi[j] + phi[m]) % 129) != ((zp[j] + 43) % 129)) lawfail++;
    }
    printf("[69] clean sections with clean reflected partner: %d (claim 1589), nontrivial %d (claim 1588)\n",
           paired, nontrivial);
    printf("[69] z(P_-j)=z(P_j): %s (%d fail)   phi_j+phi_-j = z(P_j)+43 mod 129: %s (%d fail)\n",
           zfail ? "FAIL" : "OK", zfail, lawfail ? "FAIL" : "OK", lawfail);

    /* the specific j = 780 numbers */
    if (clean[780])
        printf("[69] j=780: phi=%d (claim 23)  z(P)=%d (claim 26)  t=%ld (claim 366482)\n",
               phi[780], zp[780], tk[780L * 112]);
    if (clean[1614])
        printf("[69] j=1614: phi=%d (claim 46)  z(P)=%d\n", phi[1614], zp[1614]);

    /* how many clean j have z(P_j) = 26, i.e. how often the Theta route could apply */
    int n26 = 0;
    for (int j = 0; j < nsecs; j++) if (clean[j] && zp[j] == 26) n26++;
    printf("[69] clean sections with z(P_j)=26: %d\n", n26);

    /* dump the section table for the Python follow-up */
    FILE *f = fopen("sections69.txt", "w");
    for (int j = 0; j < nsecs; j++) fprintf(f, "%d %d %d %d\n", j, clean[j], clean[j] ? phi[j] : -1, clean[j] ? zp[j] : -1);
    fclose(f);
}

/* ---------- the terminal cone ---------- */
static int ringstate[512][MAXL];

static void terminal(int L, int *e_out, int *nout, int Umax, int Ustore[][MAXL]) {
    int x[MAXL], x0[MAXL];
    for (int i = 0; i < L; i++) x[i] = x0[i] = 1;
    long k = 0;
    for (;;) {
        Astep(x, L);
        if (x[0] == 1) { memcpy(ringstate[k % 512], x, L * sizeof(int)); k++; }
        if (!memcmp(x, x0, L * sizeof(int))) break;
    }
    long N = k;
    int n = 0;
    for (int r = 1; r <= Umax && r < N; r++) {
        int *st = ringstate[(N - 1 - r) % 512];
        /* find e with prefix sum e over first e coords and suffix = Z_{L-e} */
        int found = -1, s = 0;
        for (int e = 0; e <= L - 2; e++) {
            if (e > 0) s += st[e - 1];
            if (s != e) continue;
            int m = L - e;
            if (st[e] != 1 || st[e + 1] != m - 1) continue;
            int okz = 1;
            for (int i = e + 2; i < L; i++) if (st[i]) { okz = 0; break; }
            if (okz) { found = e; break; }
        }
        if (found < 0) break;
        e_out[n] = found;
        if (Ustore) memcpy(Ustore[n], st, found * sizeof(int));
        n++;
    }
    *nout = n;
}

int main(void) {
    build_g18();
    build_cyc13();
    run69();

    printf("\n--- terminal cone ---\n");
    static int Uref[160][MAXL];
    int eref[160], nref = 0;
    int Ls[] = {25, 29, 31, 35, 37, 39, 41, 43, 69};
    for (unsigned li = 0; li < sizeof(Ls) / sizeof(Ls[0]); li++) {
        int L = Ls[li];
        static int U[160][MAXL];
        int e[160], n;
        terminal(L, e, &n, 130, U);
        /* increments against the L=18 gap word, up to the right boundary:
           the law is only claimed while the cut has not reached the end. */
        int pred = 0, R = 0;
        while (R + 1 < 160 && pred + g18[R % ng18] <= L - 2) { pred += g18[R % ng18]; R++; }
        int incok = 1, firstfail = -1;
        for (int r = 1; r < n; r++)
            if (e[r] - e[r - 1] != g18[(r - 1) % ng18]) { incok = 0; firstfail = r + 1; break; }
        printf("       [depth %d, law predicted valid to r=%d, first deviation at r=%d]\n",
               n, R + 1, firstfail);
        int ufail = -1;
        if (nref == 0) { memcpy(Uref, U, sizeof U); memcpy(eref, e, sizeof e); nref = n; }
        else {
            int lim = n < nref ? n : nref;
            for (int r = 0; r < lim; r++) {
                if (e[r] != eref[r]) { ufail = r; break; }
                if (memcmp(U[r], Uref[r], e[r] * sizeof(int))) { ufail = r; break; }
            }
        }
        printf("L=%3d  cone depth r=1..%3d  e_1=%d  e_r increments = g18: %s   "
               "U_r matches L=25: %s\n",
               L, n, n ? e[0] : -1, incok ? "OK" : "FAIL",
               nref == n && ufail < 0 ? "OK(all)" : (ufail < 0 ? "OK(prefix)" : "FAIL"));
        if (L == 69) {
            printf("       e_r, r=1..18: ");
            for (int r = 0; r < n && r < 18; r++) printf("%d ", e[r]);
            printf("\n       r=16 -> e=%d, L-e=%d  (claim e_16=54, Z_15)\n", e[15], L - e[15]);
        }
    }
    return 0;
}
