#include <stdio.h>
#include <string.h>

typedef unsigned __int128 u128;

static inline void Astep(int *x, int L) {
    for (int i = 0; i < L - 1; i++) {
        int a = x[i];
        if (!a) continue;
        if (a == 1) { x[i] = 1 + x[i + 1]; x[i + 1] = 0; }
        else { x[i] = a - 1; x[i + 1]++; }
    }
}

static u128 basis[130];
static int nb;
static void add(u128 v) {
    for (int i = 127; i >= 0; i--) {
        if (!((v >> i) & 1)) continue;
        if (!basis[i]) { basis[i] = v; nb++; return; }
        v ^= basis[i];
    }
}

int main(void) {
    const int L = 69;
    int x[69], x0[69];
    for (int i = 0; i < L; i++) x[i] = x0[i] = 1;

    /* --- affine rank of parity vectors over the whole H_68 orbit --- */
    u128 v0 = 0;
    for (int i = 1; i < L; i++) v0 |= (u128)(x[i] & 1) << (i - 1);
    memset(basis, 0, sizeof basis); nb = 0;
    long ncp = 0;
    for (;;) {
        Astep(x, L);
        if (x[0] == 1) {
            u128 v = 0;
            for (int i = 1; i < L; i++) v |= (u128)(x[i] & 1) << (i - 1);
            add(v ^ v0);
            ncp++;
        }
        if (!memcmp(x, x0, sizeof x)) break;
    }
    printf("H_68 orbit: %ld checkpoints, affine dim of parity vectors in F_2^68 = %d "
           "(=> %d independent affine constraints)\n", ncp, nb, 68 - nb);

    /* --- same at the 16-clock section, on the 53-coordinate suffix --- */
    for (int i = 0; i < L; i++) x[i] = 1;
    memset(basis, 0, sizeof basis); nb = 0;
    int r = 53;
    u128 w0 = 0;
    for (int i = 0; i < r; i++) w0 |= (u128)(x[16 + i] & 1) << i;
    int y0[53]; memcpy(y0, x + 16, sizeof y0);
    long nouter = 0;
    int pair_ok = 1, y01_ok = 1, tail_ok = 1;
    for (;;) {
        u128 w = 0;
        for (int i = 0; i < r; i++) w |= (u128)(x[16 + i] & 1) << i;
        add(w ^ w0);
        nouter++;
        for (int j = 0; j < 20; j++)
            if (((x[16 + 2 * j] ^ x[16 + 2 * j + 1]) & 1)) pair_ok = 0;
        if (!(x[16] & 1) || !(x[17] & 1)) y01_ok = 0;
        int t = 0;
        for (int i = 40; i < 53; i++) t ^= (x[16 + i] & 1);
        if (!t) tail_ok = 0;

        int cp = 0;
        while (cp < 112) { Astep(x, L); if (x[0] == 1) cp++; }
        if (!memcmp(x + 16, y0, sizeof y0)) break;
    }
    printf("16-clock section: %ld outer states, affine dim of suffix parity in F_2^53 = %d "
           "(=> %d constraints)\n", nouter, nb, 53 - nb);
    printf("  y_{2j}=y_{2j+1} mod 2 for j<20 : %d\n", pair_ok);
    printf("  y_0=y_1=1 mod 2               : %d\n", y01_ok);
    printf("  tail y_40..y_52 odd parity    : %d\n", tail_ok);
    return 0;
}
