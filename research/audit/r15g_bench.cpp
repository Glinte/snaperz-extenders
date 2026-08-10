// UNCOMMITTED scratch benchmark -- does the R_r renormalization pay?
//
// Build:
//   g++ -O3 -march=native -std=c++17 -o /tmp/r15g research/audit/r15g_bench.cpp
//
// Three questions:
//   1. Does the renormalized scalar counter agree with blockless::run_scalar
//      on B_L, eps_L, p_L and T_L?  (eps matters: without it there is no T_L.)
//   2. How much does it actually save on the scalar path, against the predicted
//      (L-1)/(L-17)?
//   3. Is the AVX2 wavefront's cost really independent of L?  If it is, the
//      renormalization cannot help it, because the sweep count is unchanged.
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>

#include "../../src/blockless.h"

using snaperz::blockless::Result;
using Clock = std::chrono::steady_clock;

static double secs_since(Clock::time_point t0)
{
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

// ---------------------------------------------------------------- renormalized
// One macro on the r-cell tail, r = L - 16:
//     435 sweeps                      (the tail is autonomous here)
//     a_0 := 0, one sweep             (delete root)
//     30 sweeps                       (forest clock)
//     c times: a_0 += 1, one sweep    (reinject, c = a_0 after the 435)
// which is 466 + c sweeps of the full L-cell system, and lands back on a
// shielded state.  Iterate from E_r until E_r recurs: that is p_L.
//
// `skip_tail` mirrors the trailing-zero skip in blockless::run_scalar so the
// comparison is like for like.
template<bool kSkipTail>
static Result renorm(uint32_t L)
{
  const uint32_t r = L - 16;
  std::vector<uint32_t> a(r, 1);
  const uint32_t top = r - 1;
  uint32_t last = top;
  uint64_t sweeps = 0;
  int eps = 0;

  auto one_sweep = [&]() {
    if (a[top]) eps ^= 1;
    const uint32_t stop = kSkipTail ? std::min(last + 1, top) : top;
    for (uint32_t k = 0; k < stop; k++)
    {
      const uint32_t c = a[k];
      if (c == 0) continue;
      if (c == 1) { a[k] = 1 + a[k + 1]; a[k + 1] = 0; }
      else        { a[k] = c - 1;        a[k + 1] += 1; }
    }
    sweeps++;
    if (kSkipTail)
    {
      if (last < top && a[last + 1]) last++;
      while (last > 0 && a[last] == 0) last--;
    }
  };

  for (;;)
  {
    for (uint32_t i = 0; i < 435; i++) one_sweep();
    const uint32_t c = a[0];
    a[0] = 0;
    if (kSkipTail) { last = top; while (last > 0 && a[last] == 0) last--; }
    one_sweep();
    for (uint32_t i = 0; i < 30; i++) one_sweep();
    for (uint32_t i = 0; i < c; i++)
    {
      a[0] += 1;
      if (kSkipTail && last == 0 && a[0]) { /* last stays 0 */ }
      one_sweep();
    }
    bool at_E = true;
    for (uint32_t i = 0; i < r; i++) if (a[i] != 1) { at_E = false; break; }
    if (at_E) break;
  }

  Result res;
  res.p = sweeps;                       // one full cycle E_L -> E_L
  res.sweeps = res.p - (L - 1);         // B_L
  res.eps = eps;
  res.total = static_cast<uint64_t>(1 + res.eps) * res.p - (L - 1);
  return res;
}

// ------------------------------------------------------------------- avx2 probe
template<uint32_t L>
static void avx2_probe()
{
  if constexpr (snaperz::blockless::Config<L>::kFits)
  {
    const auto t0 = Clock::now();
    const Result r = snaperz::blockless::run_avx2<L>();
    const double s = secs_since(t0);
    printf("  L=%3u  fits=yes  B_L=%12llu  %8.3fs  %7.2f ns/sweep\n",
           L, (unsigned long long)r.sweeps, s, 1e9 * s / (double)r.sweeps);
  }
  else
  {
    printf("  L=%3u  fits=NO   -> falls back to the scalar path\n", L);
  }
}

int main(int argc, char **argv)
{
  // r15g renorm <L>  -- just run the renormalized counter, for cross-checking
  // a single length against another implementation.
  if (argc >= 3 && std::string(argv[1]) == "renorm")
  {
    const uint32_t L = (uint32_t)atoi(argv[2]);
    const auto t0 = Clock::now();
    const Result n = renorm<true>(L);
    printf("renorm L=%u  B_L=%llu  eps=%d  p_L=%llu  T_L=%llu  in %.1fs\n", L,
           (unsigned long long)n.sweeps, n.eps, (unsigned long long)n.p,
           (unsigned long long)n.total, secs_since(t0));
    return 0;
  }

  const bool quick = (argc > 1 && std::string(argv[1]) == "quick");

  printf("=== 1. agreement with blockless::run_scalar ===\n");
  printf("   L      B_L direct   B_L renorm   eps d/r   T_L match\n");
  for (uint32_t L : {17u, 18u, 19u, 20u, 21u, 22u, 23u, 24u, 25u, 26u, 27u,
                     28u, 29u, 30u, 31u, 32u, 33u, 34u, 35u, 36u, 40u, 44u})
  {
    const Result d = snaperz::blockless::run_scalar(L);
    const Result n = renorm<true>(L);
    const Result m = renorm<false>(L);
    printf("  %2u  %12llu %12llu    %d/%d     %s%s\n", L,
           (unsigned long long)d.sweeps, (unsigned long long)n.sweeps,
           d.eps, n.eps,
           (d.sweeps == n.sweeps && d.eps == n.eps && d.total == n.total) ? "OK" : "MISMATCH",
           (m.sweeps == n.sweeps && m.eps == n.eps) ? "" : "  (skip-tail differs!)");
  }

  printf("\n=== 2. scalar timing, direct vs renormalized ===\n");
  printf("   L      B_L        direct     renorm    speedup   predicted (L-1)/(L-17)\n");
  for (uint32_t L : {31u, 35u, 41u, 43u, 45u, 49u, 60u, 64u, 66u, 69u})
  {
    if (quick && L > 49) continue;
    auto t0 = Clock::now();
    const Result d = snaperz::blockless::run_scalar(L);
    const double td = secs_since(t0);
    t0 = Clock::now();
    const Result n = renorm<true>(L);
    const double tn = secs_since(t0);
    printf("  %2u  %12llu  %8.3fs  %8.3fs   %6.2fx    %6.2fx   %s\n", L,
           (unsigned long long)d.sweeps, td, tn, td / tn,
           (double)(L - 1) / (double)(L - 17),
           (d.sweeps == n.sweeps && d.eps == n.eps) ? "" : "MISMATCH");
    fflush(stdout);
  }

  printf("\n=== 3. is the AVX2 wavefront's cost independent of L? ===\n");
  avx2_probe<21>();
  avx2_probe<31>();
  avx2_probe<41>();
  avx2_probe<49>();
  avx2_probe<53>();
  avx2_probe<60>();
  avx2_probe<63>();
  avx2_probe<64>();
  avx2_probe<65>();
  avx2_probe<66>();
  avx2_probe<69>();
  return 0;
}
