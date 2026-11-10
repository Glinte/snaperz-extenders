// Deterministic arbitrary-state checks against sequential gates. Checks stay
// active under NDEBUG and include pipeline/lane boundaries and scalar fallback.
#include "../../src/blockless_half.h"
#include <cstdio>
#include <random>
#include <stdexcept>
using namespace snaperz::blockless::half;

static void require(bool condition, const char* message)
{
  if (!condition) throw std::runtime_error(message);
}

static void powers_and_counts()
{
  std::mt19937 rng(123);
  for (unsigned n = 1; n <= 130; ++n)
    for (unsigned count : {0, 1, 2, 15, 31, 435})
      for (unsigned start : {0, 1})
      {
        State a(n);
        for (unsigned i = 0; i < std::min(n, 90u); ++i) ++a[rng() % n];
        a[0] += n / 3;
        for (bool backward : {false, true})
        {
          auto x = a, y = a;
          for (unsigned i = 0; i < count; ++i)
            if (backward) inverse(x, start); else sweep(x, start);
          if (backward) power_simd<true>(y, count, start);
          else power_simd(y, count, start);
          if (x != y) {
            std::fprintf(stderr, "power mismatch n=%u count=%u start=%u inverse=%d\n", n,count,start,backward);
            for (unsigned j=0;j<n;++j) if(x[j]!=y[j]) std::fprintf(stderr,"  [%u] expected %u actual %u\n",j,x[j],y[j]);
          }
          require(x == y, "finite sweep block mismatch");
        }
      }
  for (unsigned n = 17; n <= 44; ++n)
  {
    auto x = run<false>(n), y = run(n);
    require(x.result.total == y.result.total && x.result.p == y.result.p &&
            x.result.eps == y.result.eps && x.result.sweeps == y.result.sweeps &&
            x.steps == y.steps && x.rewind == y.rewind && x.edge == y.edge,
            "half-cycle certificate mismatch");
  }
  puts("PASS forward/inverse powers and complete scalar/SIMD certificates 17..44");
}

static void reflection_and_injection()
{
  std::mt19937 rng(321);
  for (unsigned n = 1; n <= 129; ++n)
    for (unsigned mass = 0; mass <= 128; ++mass)
    {
      State a(n + 1), expected(n + 1), actual(n + 1);
      a[0] = expected[0] = 9;
      for (unsigned i = 0; i < mass; ++i) ++a[1 + rng() % n];
      auto x = a;
      for (unsigned rank = 1; rank <= mass; ++rank)
      {
        unsigned total = 0;
        for (unsigned i = 1; i <= n; ++i)
        {
          total += x[i];
          if (total >= rank) { ++expected[i]; break; }
        }
        if (rank < mass) inverse(x, 1);
      }
      const bool supported = retiming_simd(a, actual);
      require(!supported || expected == actual, "particle retiming mismatch");
#if SNAPERZ_BLOCKLESS_HAVE_AVX512
      require(supported == (n <= 128 && mass <= 127), "AVX512 reflection bounds");
#else
      require(supported == (n <= 64 && mass <= 127), "AVX2 reflection bounds");
#endif
      for (unsigned count : {0, 1, 3, 31})
      {
        auto u = a, v = a;
        for (unsigned i = 0; i < count; ++i) { ++u[0]; sweep(u); }
        power_simd<false, true>(v, count);
        require(u == v, "injection mismatch");
      }
    }
  puts("PASS 16,641 retiming cases and 66,564 injection cases, including fallbacks");
}

int main()
{
  powers_and_counts();
  reflection_and_injection();
}
