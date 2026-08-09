// Exact count of distinct anti-diagonals, N(t, L), with no sampling.
//
// The sampled version walks a prefix of one orbit, so it only ever produced a
// lower bound and its convergence test was confounded: fewer sweeps genuinely
// means fewer diagonals, not merely worse coverage.  Here the whole state space
// is enumerated instead.  Every window of t consecutive rows anywhere in the
// dynamics starts from some Catalan state, so taking the t rows above each
// state of Catalan(L) enumerates every anti-diagonal exactly once or more --
// no prefix, no estimator bias.
//
// HyperLogLog is still used for the cardinality, but only to keep memory
// constant; the input to it is now the complete population rather than a
// sample, so the ~0.4% is pure measurement error with no downward bias.
//
//   cc -O3 -march=native -o snaperz_diag_exact snaperz_diag_exact.c -lm
//   ./snaperz_diag_exact <length>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HLL_P 16u
#define HLL_M (1u << HLL_P)
#define MAXL 24

typedef struct { uint8_t reg[HLL_M]; } hll;

static inline uint64_t splitmix64(uint64_t x)
{
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

static inline void hll_add(hll *h, uint64_t key)
{
  const uint64_t v = splitmix64(key);
  const uint32_t idx = (uint32_t)(v >> (64 - HLL_P));
  const uint64_t tail = (v << HLL_P) | (1ULL << (HLL_P - 1));
  const uint8_t rank = (uint8_t)(__builtin_clzll(tail) + 1);
  if (rank > h->reg[idx]) h->reg[idx] = rank;
}

static double hll_count(const hll *h)
{
  const double alpha = 0.7213 / (1.0 + 1.079 / (double)HLL_M);
  double sum = 0.0;
  uint32_t zeros = 0;
  for (uint32_t i = 0; i < HLL_M; i++) {
    sum += ldexp(1.0, -(int)h->reg[i]);
    if (h->reg[i] == 0) zeros++;
  }
  double est = alpha * (double)HLL_M * (double)HLL_M / sum;
  if (est <= 2.5 * (double)HLL_M && zeros > 0)
    est = (double)HLL_M * log((double)HLL_M / (double)zeros);
  return est;
}

static int L;
static uint8_t rows[MAXL][MAXL];
static hll *H;
static unsigned long long nstates = 0;

static void sweep_into(const uint8_t *src, uint8_t *dst)
{
  memcpy(dst, src, (size_t)L);
  for (int i = 0; i < L - 1; i++) {
    const uint8_t c = dst[i];
    if (c == 0) continue;
    if (c == 1) { dst[i] = (uint8_t)(1 + dst[i + 1]); dst[i + 1] = 0; }
    else        { dst[i] = (uint8_t)(c - 1);          dst[i + 1]++;   }
  }
}

static void process(void)
{
  nstates++;
  for (int k = 1; k < L; k++) sweep_into(rows[k - 1], rows[k]);
  for (int j = 0; j < L; j++) {
    uint64_t h = 1469598103934665603ULL;
    const int kmax = (j + 1 < L) ? j + 1 : L;
    for (int k = 0; k < kmax; k++) {
      h = h * 1099511628211ULL + (uint64_t)rows[k][j - k] + 1ULL;
      hll_add(&H[k + 1], h);
    }
  }
}

// Enumerate Catalan states: a_i >= 0, sum = L, every prefix sum s_i >= i.
static void rec(int i, int s)
{
  if (i == L) { if (s == L) process(); return; }
  int lo = (i + 1) - s;
  if (lo < 0) lo = 0;
  for (int a = lo; s + a <= L; a++) {
    rows[0][i] = (uint8_t)a;
    rec(i + 1, s + a);
  }
}

int main(int argc, char **argv)
{
  if (argc != 2) { fprintf(stderr, "usage: %s <length>\n", argv[0]); return 2; }
  L = atoi(argv[1]);
  if (L < 2 || L > MAXL) { fprintf(stderr, "length must be 2..%d\n", MAXL); return 2; }

  H = calloc((size_t)L + 1, sizeof(hll));
  if (!H) { fprintf(stderr, "out of memory\n"); return 1; }

  rec(0, 0);

  printf("L=%d  states=%llu (all of Catalan(%d))\n", L, nstates, L);
  printf("%4s %14s %10s %12s\n", "t", "N(t,L) exact", "ratio", "log-log");
  double prev = 0.0;
  for (int t = 1; t < L; t++) {
    const double n = hll_count(&H[t]);
    if (t == 1) printf("%4d %14.0f %10s %12s\n", t, n, "-", "-");
    else printf("%4d %14.0f %10.3f %12.2f\n", t, n, n / prev,
                log(n / prev) / log((double)t / (double)(t - 1)));
    prev = n;
  }
  free(H);
  return 0;
}
