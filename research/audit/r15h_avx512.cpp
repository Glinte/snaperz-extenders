// Round-15 audit, part H: validates the 512-bit blockless windows in
// src/blockless_avx512.h against the 256-bit ones and against the B_L table.
//
//   g++ -O3 -march=native -std=c++17 -o build/r15h audit/r15h_avx512.cpp
//
// blockless.h's 256-bit wavefront holds two 32-lane uint8 windows, so it fits
// only kSegCount <= 64, i.e. L <= 64; past that the counter dropped to the
// scalar sweep. That cliff, not the cost per cell, is what made lengths just
// past 64 expensive -- the wavefront costs ~10ns per sweep independently of L,
// which is also why the R_r renormalization cannot help it (r15a..r15g): R_r
// removes cells, and the wavefront does not care how many cells there are.
//
// Widening the windows to 64 lanes reaches L <= 127 instead. This checks it.
#include <cstdio>
#include <cstdint>
#include <chrono>
#include <string>

#include "../../src/blockless.h"

using Clock = std::chrono::steady_clock;
using snaperz::blockless::Result;

static double secs_since(Clock::time_point t0)
{
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

#if !SNAPERZ_BLOCKLESS_HAVE_AVX512
int main(void)
{
  printf("built without AVX512BW; nothing to check\n");
  return 0;
}
#else

// B_L for L = 4..46, from docs/research-notes.md.
static const uint64_t kB[] = {
  0, 0, 0, 0, 3, 9, 9, 22, 23, 23, 23, 54, 23, 117, 185, 211, 451, 451, 451,
  918, 451, 1854, 919, 1853, 919, 3729, 4667, 5608, 7493, 11251, 4673, 115288,
  6559, 36504, 8887, 201627, 37551, 386120, 15927, 705934, 21567, 382512,
  185403, 2232272, 38513, 729845, 33767 };

// Every L where both widths work must agree on all four outputs.
template<uint32_t L>
struct CrossCheck
{
  static void run(int &bad, int &n)
  {
    CrossCheck<L - 1>::run(bad, n);
    if constexpr (snaperz::blockless::Config<L>::kFits &&
                  snaperz::blockless::wide::Config<L>::kFits)
    {
      const Result a = snaperz::blockless::run_avx2<L>();
      const Result b = snaperz::blockless::wide::run<L>();
      n++;
      const bool agree = a.sweeps == b.sweeps && a.eps == b.eps &&
                         a.p == b.p && a.total == b.total;
      const bool tabled = (L < sizeof(kB) / sizeof(kB[0])) && kB[L];
      if (!agree || (tabled && a.sweeps != kB[L]))
      {
        bad++;
        printf("  L=%3u  256-bit B=%llu eps=%d   512-bit B=%llu eps=%d   %s\n", L,
               (unsigned long long)a.sweeps, a.eps,
               (unsigned long long)b.sweeps, b.eps,
               agree ? "disagrees with the notes" : "MISMATCH");
      }
    }
  }
};
template<> struct CrossCheck<16> { static void run(int &, int &) {} };

template<uint32_t L>
static void report(void)
{
  if constexpr (!snaperz::blockless::wide::Config<L>::kFits)
  {
    printf("  L=%3u  past the 512-bit limit too; scalar path\n", L);
  }
  else
  {
    const auto t0 = Clock::now();
    const Result r = snaperz::blockless::wide::run<L>();
    const double s = secs_since(t0);
    const char *verdict = "";
    if (L < sizeof(kB) / sizeof(kB[0]) && kB[L])
      verdict = (r.sweeps == kB[L]) ? "matches the notes" : "*** WRONG ***";
    printf("  L=%3u  B_L=%12llu  eps=%d  %8.3fs  %6.2f ns/sweep  %s\n", L,
           (unsigned long long)r.sweeps, r.eps, s,
           r.sweeps ? 1e9 * s / (double)r.sweeps : 0.0, verdict);
  }
  fflush(stdout);
}

int main(void)
{
  printf("--- 256-bit vs 512-bit windows, every L in 17..46 ---\n");
  int bad = 0, n = 0;
  CrossCheck<46>::run(bad, n);
  printf("  %d lengths compared, %d mismatches\n\n", n, bad);

  printf("--- the band the 256-bit windows cannot reach ---\n");
  report<60>();   // still narrow, for scale
  report<64>();   // last length the narrow windows hold
  report<66>();
  report<69>();

  printf("\n--- narrow is faster where it applies, so it stays first ---\n");
  {
    auto t0 = Clock::now();
    const Result a = snaperz::blockless::run_avx2<64>();
    const double ta = secs_since(t0);
    t0 = Clock::now();
    const Result b = snaperz::blockless::wide::run<64>();
    const double tb = secs_since(t0);
    printf("  L=64  256-bit %.3fs   512-bit %.3fs   %.2fx   agree=%s\n",
           ta, tb, ta / tb, (a.sweeps == b.sweeps) ? "yes" : "NO");
  }
  return bad ? 1 : 0;
}
#endif
