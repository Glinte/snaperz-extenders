// Counts distinct anti-diagonals of the Snaperz space-time diagram.
//
// To evaluate A^t in one left-to-right pass you must carry the anti-diagonal
// (a^(1)_j, a^(2)_{j-1}, ..., a^(t)_{j-t+1}), so the number of distinct
// anti-diagonals lower-bounds the size of any transducer for A^t.  If that
// count is polynomial in t then A^t has a compact representation and repeated
// squaring gives B_L in polylog time; if it is exponential, every fast-forward
// route that needs to compose A with itself is closed.
//
// Python could only reach ~4.5e6 samples, which was far short of saturation and
// bent the growth curve downward three separate ways.  This reaches ~1e9 using
// HyperLogLog for the cardinalities, so the memory cost is constant in the
// number of distinct diagonals rather than linear.
//
//   cc -O3 -march=native -o snaperz_diag snaperz_diag.c -lm
//   ./snaperz_diag <length> <sweeps> <max_t>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HLL_P 16u
#define HLL_M (1u << HLL_P)

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
  // Guarantee a set bit so the rank is bounded even when the tail is all zero.
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
  // Small-range correction: below ~2.5m the raw estimator is badly biased.
  if (est <= 2.5 * (double)HLL_M && zeros > 0)
    est = (double)HLL_M * log((double)HLL_M / (double)zeros);
  return est;
}

// One ordinary blockless sweep: (0,y)->(0,y), (1,y)->(1+y,0), (x,y)->(x-1,y+1).
static void sweep(uint32_t *a, int len)
{
  for (int i = 0; i < len - 1; i++) {
    const uint32_t c = a[i];
    if (c == 0) continue;
    if (c == 1) { a[i] = 1 + a[i + 1]; a[i + 1] = 0; }
    else        { a[i] = c - 1;        a[i + 1] += 1; }
  }
}

int main(int argc, char **argv)
{
  if (argc != 4) { fprintf(stderr, "usage: %s <length> <sweeps> <max_t>\n", argv[0]); return 2; }
  const int len = atoi(argv[1]);
  const long sweeps = atol(argv[2]);
  const int max_t = atoi(argv[3]);
  if (len < max_t + 2 || max_t < 1) { fprintf(stderr, "bad parameters\n"); return 2; }

  const int rows = max_t + 1;
  uint32_t *ring = calloc((size_t)rows * len, sizeof(uint32_t));
  uint32_t *cur = malloc((size_t)len * sizeof(uint32_t));
  hll *H = calloc((size_t)max_t + 1, sizeof(hll));
  if (!ring || !cur || !H) { fprintf(stderr, "out of memory\n"); return 1; }

  for (int i = 0; i < len; i++) cur[i] = 1;

  unsigned long long samples = 0;
  for (long r = 0; r <= sweeps; r++) {
    memcpy(ring + (size_t)(r % rows) * len, cur, (size_t)len * sizeof(uint32_t));
    if (r >= max_t) {
      for (int j = max_t - 1; j < len; j++) {
        uint64_t h = 1469598103934665603ULL;
        for (int k = 0; k < max_t; k++) {
          const uint32_t v = ring[(size_t)((r - k) % rows) * len + (j - k)];
          h = h * 1099511628211ULL + (uint64_t)v + 1ULL;
          hll_add(&H[k + 1], h);
        }
        samples++;
      }
    }
    sweep(cur, len);
  }

  printf("L=%d sweeps=%ld samples=%llu (per t)\n", len, sweeps, samples);
  printf("%4s %16s %10s\n", "t", "distinct (HLL)", "per-level");
  double prev = 0.0;
  for (int t = 1; t <= max_t; t++) {
    const double n = hll_count(&H[t]);
    if (t == 1) printf("%4d %16.0f %10s\n", t, n, "-");
    else        printf("%4d %16.0f %10.3f\n", t, n, n / prev);
    prev = n;
  }
  free(ring); free(cur); free(H);
  return 0;
}
