// Correctness checks remain active in release builds. Timings are end-to-end.
#include "../../src/blockless_half.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
using namespace snaperz::blockless;
using Clock = std::chrono::steady_clock;

static void equal(const Result& a, const Result& b)
{
  if (a.sweeps != b.sweeps || a.eps != b.eps || a.p != b.p || a.total != b.total)
    throw std::runtime_error("counter mismatch");
}

template<unsigned L> void bench()
{
  auto start = Clock::now();
  auto h = half::run(L);
  double th = std::chrono::duration<double>(Clock::now()-start).count();
  start = Clock::now();
  auto s = run_scalar(L);
  double ts = std::chrono::duration<double>(Clock::now()-start).count();
  equal(h.result, s);
  start = Clock::now();
  Result v;
  if constexpr (Config<L>::kFits) v = run_avx2<L>();
#if SNAPERZ_BLOCKLESS_HAVE_AVX512
  else v = wide::run<L>();
#else
  else v = s; // No SIMD implementation available at this width.
#endif
  double tv = std::chrono::duration<double>(Clock::now()-start).count();
#if !SNAPERZ_BLOCKLESS_HAVE_AVX512
  if constexpr (!Config<L>::kFits) tv = -1; // Unavailable, not a timing.
#endif
  equal(h.result, v);
  printf("L=%u pulses=%llu steps=%llu %s rewind=%llu half=%.6fs scalar=%.6fs simd=%.6fs scalar/half=%.2f simd/half=%.2f\n",
         L, (unsigned long long)h.result.total, (unsigned long long)h.steps,
         h.edge ? "edge" : "vertex", (unsigned long long)h.rewind, th, ts, tv, ts/th, tv/th);
  fflush(stdout);
}

int main(int argc, char** argv)
{
  if (argc > 1 && std::string(argv[1]) == "66") { bench<66>(); return 0; }
  bool edge = false, vertex = false, odd = false, eps = false;
  for (unsigned n = 17; n <= 50; ++n)
  {
    auto h = half::run(n);
    equal(h.result, run_scalar(n));
    printf("verified L=%u\n", n); fflush(stdout);
    edge |= h.edge; vertex |= !h.edge; odd |= h.result.p & 1; eps |= h.result.eps;
  }
  if (!(edge && vertex && odd && eps)) throw std::runtime_error("missing branch coverage");
  auto h = half::run(50);
  if (h.steps != 468 || h.rewind != 15 || h.result.p != 439822 || h.result.total != 439773)
    throw std::runtime_error("length-50 certificate mismatch");
  puts("PASS: L=17..50, both meeting types, odd/even period, nonzero epsilon, L=50 certificate");
  bench<49>(); bench<50>(); bench<60>(); bench<64>();
  if (argc > 1) { bench<65>(); bench<66>(); bench<69>(); }
}
