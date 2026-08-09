/* Branchless + vectorised variants.  The L=49 loop is branch-mispredict bound,
   which is the thing the cumulative representation is actually able to fix. */
#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct { uint64_t q, p; int parity; } Res;

/* --- B: a-coordinates, branchy, carried register (best from round 1) ----- */
static Res run_b(int L) {
    uint8_t *s = calloc((size_t)L + 2, 1);
    for (int i = 0; i < L; i++) s[i] = 1;
    const int top = L - 1;
    int last = top, parity = 0; uint64_t b = 0, q = 0;
    for (;;) {
        if (last == top) parity ^= 1;
        if (s[0] == 1) q++;
        int stop = last + 1; if (stop > top) stop = top;
        int c = s[0];
        for (int i = 0; i < stop; i++) {
            int n = s[i + 1];
            if (c == 0)      { s[i] = 0;                c = n; }
            else if (c == 1) { s[i] = (uint8_t)(1 + n); c = 0; }
            else             { s[i] = (uint8_t)(c - 1); c = n + 1; }
        }
        s[stop] = (uint8_t)c;
        b++;
        if (last < top && s[last + 1]) last++;
        while (!s[last]) last--;
        if (s[0] == L) break;
    }
    free(s);
    Res r = { q, b + (uint64_t)(L - 1), parity }; return r;
}

/* --- F: a-coordinates, branchless (cmov) -------------------------------- */
static Res run_f(int L) {
    uint8_t *s = calloc((size_t)L + 2, 1);
    for (int i = 0; i < L; i++) s[i] = 1;
    const int top = L - 1;
    int last = top, parity = 0; uint64_t b = 0, q = 0;
    for (;;) {
        if (last == top) parity ^= 1;
        if (s[0] == 1) q++;
        int stop = last + 1; if (stop > top) stop = top;
        int c = s[0];
        for (int i = 0; i < stop; i++) {
            int n = s[i + 1];
            int w  = (c == 0) ? 0 : ((c == 1) ? 1 + n : c - 1);
            int cn = (c == 0) ? n : ((c == 1) ? 0     : n + 1);
            s[i] = (uint8_t)w;
            c = cn;
        }
        s[stop] = (uint8_t)c;
        b++;
        if (last < top && s[last + 1]) last++;
        while (!s[last]) last--;
        if (s[0] == L) break;
    }
    free(s);
    Res r = { q, b + (uint64_t)(L - 1), parity }; return r;
}

#define SM1(i) ((i) ? s[(i) - 1] : 0)

/* --- E: cumulative, branchless (cmov) ----------------------------------- */
static Res run_e(int L) {
    uint8_t *s = calloc((size_t)L + 2, 1);
    for (int i = 0; i < L; i++) s[i] = (uint8_t)(i + 1);
    const int top = L - 1;
    int last = top, parity = 0; uint64_t b = 0, q = 0;
    for (;;) {
        if (last == top) parity ^= 1;
        if (s[0] == 1) q++;
        int stop = last + 1; if (stop > top) stop = top;
        int prev = 0;
        for (int i = 0; i < stop; i++) {
            int cur = s[i], nxt = s[i + 1];
            int a = cur - prev;
            int r = (a == 0) ? cur : ((a == 1) ? nxt : cur - 1);
            s[i] = (uint8_t)r;
            prev = r;
        }
        b++;
        if (last < top && s[last + 1] > s[last]) last++;
        while (s[last] == SM1(last)) last--;
        if (s[0] == L) break;
    }
    free(s);
    Res r = { q, b + (uint64_t)(L - 1), parity }; return r;
}

/* --- G: cumulative + 4-cell table carry scan + SSE passes 1 and 3 -------- */
/* carry state 0=Z (keep), 1=D (s_i-1), 2=W (s_i = s_{i+1}) */
static const uint8_t DELTA[3][3] = { {0,2,1}, {2,1,1}, {0,0,0} };
static uint32_t TAB[3][256];   /* [carry_in][4 packed symbols] -> 4 state bytes */

static void build_tables(void) {
    for (int st0 = 0; st0 < 3; st0++)
        for (int idx = 0; idx < 256; idx++) {
            int st = st0; uint32_t out = 0;
            for (int k = 0; k < 4; k++) {
                int sym = ((idx >> k) & 1) | (((idx >> (4 + k)) & 1) << 1);
                st = DELTA[st][sym];
                out |= (uint32_t)st << (8 * k);
            }
            TAB[st0][idx] = out;
        }
}

static Res run_g(int L) {
    int cap = ((L + 80) & ~15);
    uint8_t *raw = aligned_alloc(64, (size_t)cap + 64);
    memset(raw, 0, (size_t)cap + 64);
    uint8_t *s = raw + 16;                 /* s[-1] readable, holds 0 */
    uint8_t *sym = aligned_alloc(64, (size_t)cap + 64);
    uint8_t *st  = aligned_alloc(64, (size_t)cap + 64);
    memset(sym, 0, (size_t)cap + 64); memset(st, 0, (size_t)cap + 64);
    for (int i = 0; i < L; i++) s[i] = (uint8_t)(i + 1);
    for (int i = L; i < cap; i++) s[i] = (uint8_t)L;

    const __m128i two = _mm_set1_epi8(2), one = _mm_set1_epi8(1);
    const __m128i c1 = _mm_set1_epi8(1), c2 = _mm_set1_epi8(2), zero = _mm_setzero_si128();
    const int top = L - 1;
    int last = top, parity = 0; uint64_t b = 0, q = 0;
    for (;;) {
        if (last == top) parity ^= 1;
        if (s[0] == 1) q++;
        int stop = last + 1; if (stop > top) stop = top;
        int nb = (stop + 3) >> 2;              /* 4-cell carry blocks */
        int nsym = nb * 4, nvec = (nsym + 15) >> 4;

        /* pass 1 (SIMD): sym_i = min(s_i - s_{i-1}, 2), and pack to 4-bit masks */
        uint32_t m0[8], m1[8];
        for (int v = 0; v < nvec; v++) {
            __m128i cur = _mm_loadu_si128((const __m128i *)(s + v * 16));
            __m128i pre = _mm_loadu_si128((const __m128i *)(s + v * 16 - 1));
            __m128i a = _mm_min_epu8(_mm_sub_epi8(cur, pre), two);
            m0[v] = (uint32_t)_mm_movemask_epi8(_mm_slli_epi16(a, 7));
            m1[v] = (uint32_t)_mm_movemask_epi8(_mm_slli_epi16(a, 6));
        }
        /* pass 2 (serial, branchless): one L1 table lookup per 4 cells */
        uint32_t cs = 0;
        for (int blk = 0; blk < nb; blk++) {
            int v = blk >> 2, sh = (blk & 3) * 4;
            uint32_t idx = ((m0[v] >> sh) & 0xF) | (((m1[v] >> sh) & 0xF) << 4);
            uint32_t out = TAB[cs][idx];
            *(uint32_t *)(st + blk * 4) = out;
            cs = out >> 24;
        }
        _mm_storeu_si128((__m128i *)(st + stop), zero);   /* kill tail carries */

        /* pass 3 (SIMD, branchless): s_i' from OLD s and the carry state */
        for (int v = 0; v < nvec; v++) {
            __m128i cur = _mm_loadu_si128((const __m128i *)(s + v * 16));
            __m128i nxt = _mm_loadu_si128((const __m128i *)(s + v * 16 + 1));
            __m128i q3 = _mm_loadu_si128((const __m128i *)(st + v * 16));
            __m128i r = _mm_blendv_epi8(cur, _mm_sub_epi8(cur, one), _mm_cmpeq_epi8(q3, c1));
            r = _mm_blendv_epi8(r, nxt, _mm_cmpeq_epi8(q3, c2));
            _mm_storeu_si128((__m128i *)(s + v * 16), r);
        }
        b++;
        if (last < top && s[last + 1] > s[last]) last++;
        while (s[last] == SM1(last)) last--;
        if (s[0] == L) break;
    }
    free(raw); free(sym); free(st);
    Res r = { q, b + (uint64_t)(L - 1), parity }; return r;
}

static double now(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    build_tables();
    int L = argc > 1 ? atoi(argv[1]) : 41;
    struct { const char *n; Res (*f)(int); } v[] = {
        {"B  a-coord branchy   ", run_b},
        {"F  a-coord branchless", run_f},
        {"E  cumulative brchles", run_e},
        {"G  cumulative+FSM+SSE", run_g},
    };
    Res ref; int first = 1;
    for (unsigned i = 0; i < sizeof v / sizeof *v; i++) {
        double t0 = now(); Res r = v[i].f(L); double dt = now() - t0;
        if (first) { ref = r; first = 0; }
        printf("%s  p=%-11llu q=%-10llu par=%d  %8.3f s  %s\n", v[i].n,
               (unsigned long long)r.p, (unsigned long long)r.q, r.parity, dt,
               (r.p == ref.p && r.q == ref.q && r.parity == ref.parity) ? "ok" : "MISMATCH");
        fflush(stdout);
    }
    return 0;
}
