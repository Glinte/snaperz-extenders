#pragma once
#define SNAPERZ_BLOCKLESS_H_INCLUDED 1
// AVX2 wavefront for the *blockless* sweep A of a P=12 snaperz extender.
//
// The blockless state is a_0 .. a_{L-1} with sum L, starting at E_L = (1,...,1).
// One sweep is, for k = 0 .. L-2 (segment L-1 never acts):
//     a[k] == 0 : nothing
//     a[k] == 1 : a[k] += a[k+1]; a[k+1] = 0        (pull)
//     a[k]  > 1 : a[k] -= 1;     a[k+1] += 1        (push, kPushLimit == 1)
//
// Scanning from E_L to R_L = (L,0,...,0) takes B_L sweeps, and
//     p_L = B_L + (L - 1),   T_L = (1 + eps_L) * p_L - (L - 1)
// where eps_L is the parity of the sweeps whose starting state has a[L-1] != 0.
//
// Two differences from the real-extender wavefront make this cheaper: there is
// no running block counter (the only segment that must not act is the top one,
// which is positional), and the push limit is constant.
//
// The top and first segments are located with a one-hot marker register that is
// permuted in lockstep with the window, so it stays aligned with the segment it
// was seeded next to. Marker encoding: 2 at segment L-1, 1 at segment 0.

#include <immintrin.h>
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <type_traits>
#include <vector>

#include "constants.h"

namespace snaperz
{
namespace blockless
{
  static constexpr uint32_t kElems = 32;   // uint8 lanes per window

  template<uint32_t L>
  struct Config
  {
    static constexpr uint32_t kSegCount = (L > 2 * kElems) ? L : (L + (L & 1));
    static constexpr uint32_t kSatCount = (kSegCount < 2 * kElems) ? kSegCount : 2 * kElems;
    static constexpr bool kFits = (kSegCount == kSatCount);
    static constexpr uint32_t kLastElem =
      ((kElems - 1) < (kSatCount / 2 - 1)) ? (kElems - 1) : (kSatCount / 2 - 1);
  };

  static inline __m256i shr1(__m256i v)
  {
    __m256i hi = _mm256_permute2x128_si256(v, v, 0x81);
    return _mm256_alignr_epi8(hi, v, 1);
  }

  template<uint32_t LastElem>
  static inline __m256i lane_mask()
  {
    const __m256i idx = _mm256_setr_epi8(
      0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
      16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31);
    return _mm256_cmpeq_epi8(idx, _mm256_set1_epi8(static_cast<char>(LastElem)));
  }

  // Shift down one lane and re-insert the evicted lane 0 at kLastElem.
  template<uint32_t LastElem>
  static inline __m256i rot(__m256i v)
  {
    __m256i evicted = _mm256_broadcastb_epi8(_mm256_castsi256_si128(v));
    return _mm256_blendv_epi8(shr1(v), evicted, lane_mask<LastElem>());
  }

  // Shift down one lane and insert `value` at kLastElem.
  template<uint32_t LastElem>
  static inline __m256i shift_in(__m256i v, uint8_t value)
  {
    return _mm256_blendv_epi8(shr1(v), _mm256_set1_epi8(static_cast<char>(value)),
                              lane_mask<LastElem>());
  }

  struct Result
  {
    uint64_t sweeps;      // B_L
    int eps;              // eps_L
    uint64_t p;           // p_L
    uint64_t total;       // T_L
  };

  template<uint32_t L>
  Result run_avx2()
  {
    using C = Config<L>;
    static_assert(L >= 2, "degenerate");

    // Segment k is inserted at step k, into window (k & 1) ^ 1, and stays in
    // that window for good since the rotation never crosses windows. A window
    // is the current one on the steps of its own parity, so both markers are
    // live on statically known steps.
    constexpr uint32_t kFirstParity = 1;
    constexpr uint32_t kTopParity = ((L - 1) & 1u) ^ 1u;

    const __m256i zeros = _mm256_setzero_si256();
    const __m256i ones  = _mm256_set1_epi8(1);
    const __m256i twos  = _mm256_set1_epi8(2);
    const __m256i full  = _mm256_set1_epi8(static_cast<char>(L));

    __m256i w0 = zeros, w1 = zeros;    // segment values, by window
    __m256i m0 = zeros, m1 = zeros;    // one-hot markers, permuted alongside

    uint64_t sweeps = 0;
    int eps = 1;                       // sweep 0 starts at E_L, whose top is 1
    bool done = false;

    // --- warm-up: segments still stream in from the initial state ---
    uint32_t parity = 0;
    for (uint32_t step = 0; step < C::kSatCount && !done; step++)
    {
      __m256i curr = parity ? w1 : w0;
      __m256i mcur = parity ? m1 : m0;
      __m256i at_top = _mm256_cmpeq_epi8(mcur, twos);
      __m256i at_first = _mm256_cmpeq_epi8(mcur, ones);

      if (_mm256_movemask_epi8(_mm256_and_si256(at_first, _mm256_cmpeq_epi8(curr, full))))
      {
        done = true;
        break;
      }
      if (_mm256_movemask_epi8(at_first))
      {
        sweeps++;
      }
      if (_mm256_movemask_epi8(_mm256_andnot_si256(_mm256_cmpeq_epi8(curr, zeros), at_top)))
      {
        eps ^= 1;
      }

      const uint32_t p = step;
      const uint8_t value  = (p < L) ? 1 : 0;
      const uint8_t marker = (p == L - 1) ? 2 : ((p == 0) ? 1 : 0);
      __m256i next = shift_in<C::kLastElem>(parity ? w0 : w1, value);
      __m256i mnext = shift_in<C::kLastElem>(parity ? m0 : m1, marker);

      __m256i eq1 = _mm256_cmpeq_epi8(curr, ones);
      __m256i gt1 = _mm256_cmpgt_epi8(curr, ones);
      __m256i delta = _mm256_add_epi8(_mm256_and_si256(eq1, next), gt1);
      delta = _mm256_andnot_si256(at_top, delta);

      curr = _mm256_add_epi8(curr, delta);
      next = _mm256_sub_epi8(next, delta);
      if (parity) { w1 = curr; w0 = next; m0 = mnext; }
      else        { w0 = curr; w1 = next; m1 = mnext; }
      parity ^= 1u;
    }

    // --- steady state ---
    // The markers are one-hot and only ever move down one lane per rotation of
    // their own window, so their position is a scalar. Tracking it as an index
    // instead of a vector removes two shuffle chains per pulse; the tests
    // become a movemask and a shift, and the suppression mask is rebuilt from
    // the index. Seed the indices from the markers the warm-up already built.
    constexpr uint32_t kPeriod = C::kSatCount / 2;
    const __m256i lane_index = _mm256_setr_epi8(
      0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
      16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31);

    uint32_t top_lane = __builtin_ctz(static_cast<uint32_t>(_mm256_movemask_epi8(
      _mm256_cmpeq_epi8((kTopParity == 0) ? m0 : m1, twos))));
    uint32_t first_lane = __builtin_ctz(static_cast<uint32_t>(_mm256_movemask_epi8(
      _mm256_cmpeq_epi8(m1, ones))));

    while (!done)
    {
      // ---- parity 0: w0 is current, w1 is rotated ----
      {
        __m256i curr = w0;
        __m256i delta;
        if constexpr (kTopParity == 0)
        {
          eps ^= (~static_cast<uint32_t>(_mm256_movemask_epi8(
            _mm256_cmpeq_epi8(curr, zeros))) >> top_lane) & 1u;
        }
        __m256i next = rot<C::kLastElem>(w1);
        __m256i eq1 = _mm256_cmpeq_epi8(curr, ones);
        __m256i gt1 = _mm256_cmpgt_epi8(curr, ones);
        delta = _mm256_add_epi8(_mm256_and_si256(eq1, next), gt1);
        if constexpr (kTopParity == 0)
        {
          delta = _mm256_andnot_si256(
            _mm256_cmpeq_epi8(lane_index, _mm256_set1_epi8(static_cast<char>(top_lane))),
            delta);
        }
        w0 = _mm256_add_epi8(curr, delta);
        w1 = _mm256_sub_epi8(next, delta);
        // w1 rotated: anything living in it moved down a lane.
        first_lane = first_lane ? first_lane - 1 : kPeriod - 1;
        if constexpr (kTopParity == 1)
        {
          top_lane = top_lane ? top_lane - 1 : kPeriod - 1;
        }
      }

      // ---- parity 1: w1 is current, w0 is rotated ----
      {
        __m256i curr = w1;
        if ((static_cast<uint32_t>(_mm256_movemask_epi8(
              _mm256_cmpeq_epi8(curr, full))) >> first_lane) & 1u)
        {
          done = true;
          break;
        }
        sweeps++;
        if constexpr (kTopParity == 1)
        {
          eps ^= (~static_cast<uint32_t>(_mm256_movemask_epi8(
            _mm256_cmpeq_epi8(curr, zeros))) >> top_lane) & 1u;
        }
        __m256i next = rot<C::kLastElem>(w0);
        __m256i eq1 = _mm256_cmpeq_epi8(curr, ones);
        __m256i gt1 = _mm256_cmpgt_epi8(curr, ones);
        __m256i delta = _mm256_add_epi8(_mm256_and_si256(eq1, next), gt1);
        if constexpr (kTopParity == 1)
        {
          delta = _mm256_andnot_si256(
            _mm256_cmpeq_epi8(lane_index, _mm256_set1_epi8(static_cast<char>(top_lane))),
            delta);
        }
        w1 = _mm256_add_epi8(curr, delta);
        w0 = _mm256_sub_epi8(next, delta);
        if constexpr (kTopParity == 0)
        {
          top_lane = top_lane ? top_lane - 1 : kPeriod - 1;
        }
      }
    }

    Result r;
    r.sweeps = sweeps;
    r.eps = eps;
    r.p = sweeps + (L - 1);
    r.total = static_cast<uint64_t>(1 + r.eps) * r.p - (L - 1);
    return r;
  }

  // Scalar reference, used when the extender no longer fits inside the two
  // windows. Same sweep, written the obvious way; the trailing zeros are
  // skipped because everything past the last non-empty segment is fixed.
  inline Result run_scalar(uint32_t L)
  {
    std::vector<uint32_t> a(L, 1);
    const uint32_t top = L - 1;
    uint32_t last = top;
    uint64_t sweeps = 0;
    int eps = 0;
    for (;;)
    {
      if (last == top)
      {
        eps ^= 1;
      }
      if (a[0] == L)
      {
        break;
      }
      const uint32_t stop = std::min(last + 1, top);
      for (uint32_t k = 0; k < stop; k++)
      {
        const uint32_t c = a[k];
        if (c == 0)
        {
          continue;
        }
        if (c == 1)
        {
          a[k] = 1 + a[k + 1];
          a[k + 1] = 0;
        }
        else
        {
          a[k] = c - 1;
          a[k + 1] += 1;
        }
      }
      sweeps++;
      if (last < top && a[last + 1])
      {
        last++;
      }
      while (a[last] == 0)
      {
        last--;
      }
    }
    Result r;
    r.sweeps = sweeps;
    r.eps = eps;
    r.p = sweeps + (L - 1);
    r.total = static_cast<uint64_t>(1 + r.eps) * r.p - (L - 1);
    return r;
  }

} // namespace blockless
} // namespace snaperz

// Wider windows for the lengths the 32-lane ones cannot hold. Included here,
// after Result exists, because it is an alternative body for the same wavefront
// rather than a separate algorithm.
#include "blockless_avx512.h"

namespace snaperz
{
namespace blockless
{
  // Narrow windows first: they reach L <= 64 and are the fastest per sweep.
  // Wide windows carry L <= 127, at about 0.6x the per-sweep rate but against a
  // scalar path that is roughly 20x slower again, so the order matters.
  inline Result run()
  {
    if constexpr (Config<kLength>::kFits)
    {
      return run_avx2<kLength>();
    }
#if SNAPERZ_BLOCKLESS_HAVE_AVX512
    else if constexpr (wide::Config<kLength>::kFits)
    {
      return wide::run<kLength>();
    }
#endif
    else
    {
      return run_scalar(kLength);
    }
  }
} // namespace blockless
} // namespace snaperz
