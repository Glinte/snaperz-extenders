#pragma once
//
// 512-bit windows for the blockless wavefront. Included by blockless.h, which
// defines Result and picks between this, the 256-bit version and the scalar
// fallback; it is not meant to be included on its own.
//
// The 256-bit wavefront keeps the extender in two 32-lane uint8 windows, so it
// fits kSegCount <= 64, i.e. L <= 64. Past that blockless.h had to fall back to
// the scalar sweep, which is roughly 20x slower, and that cliff -- not the cost
// per cell -- is what made lengths just past 64 expensive. Widening the windows
// to 64 lanes moves the cliff to L <= 127.
//
// 127 rather than 128 because the lanes are compared as signed bytes and the
// terminal state is R_L = (L, 0, ..., 0), so L itself has to fit in an int8.
//
// This is *slower* per sweep than the 256-bit version (measured 0.61x at
// L = 64: the byte-granular shift below costs a cross-lane permute plus two
// shifts, where 256-bit gets it from one alignr, and the wider code clocks
// lower). It is therefore only ever used for lengths the narrow windows cannot
// hold, never as a replacement.

#ifndef SNAPERZ_BLOCKLESS_H_INCLUDED
#error "include blockless.h, not blockless_avx512.h"
#endif

#if defined(__AVX512BW__) && defined(__AVX512F__)
#define SNAPERZ_BLOCKLESS_HAVE_AVX512 1

#include <immintrin.h>
#include <cstdint>

namespace snaperz
{
namespace blockless
{
namespace wide
{
  static constexpr uint32_t kElems = 64;   // uint8 lanes per window

  template<uint32_t L>
  struct Config
  {
    static constexpr uint32_t kSegCount = (L > 2 * kElems) ? L : (L + (L & 1));
    static constexpr uint32_t kSatCount = (kSegCount < 2 * kElems) ? kSegCount : 2 * kElems;
    static constexpr bool kFits = (kSegCount == kSatCount) && (L <= 127);
    static constexpr uint32_t kLastElem =
      ((kElems - 1) < (kSatCount / 2 - 1)) ? (kElems - 1) : (kSatCount / 2 - 1);
  };

  // Shift the whole register down one byte, filling the top lane with zero.
  // A byte-granular permute would need AVX512VBMI, which is a later extension
  // than the rest of this file needs, so move the bytes within each qword and
  // stitch the qword boundaries back up by hand.
  static inline __m512i shr1(__m512i v)
  {
    const __m512i idx = _mm512_setr_epi64(1, 2, 3, 4, 5, 6, 7, 7);
    const __m512i nxt = _mm512_maskz_permutexvar_epi64(0x7F, idx, v);
    return _mm512_or_si512(_mm512_srli_epi64(v, 8), _mm512_slli_epi64(nxt, 56));
  }

  // Shift down one lane and re-insert the evicted lane 0 at kLastElem.
  template<uint32_t LastElem>
  static inline __m512i rot(__m512i v)
  {
    const int evicted = _mm_extract_epi8(_mm512_castsi512_si128(v), 0);
    return _mm512_mask_set1_epi8(shr1(v), 1ULL << LastElem, static_cast<char>(evicted));
  }

  // Shift down one lane and insert `value` at kLastElem.
  template<uint32_t LastElem>
  static inline __m512i shift_in(__m512i v, uint8_t value)
  {
    return _mm512_mask_set1_epi8(shr1(v), 1ULL << LastElem, static_cast<char>(value));
  }

  // Same wavefront as blockless::run_avx2, lane for lane; see the commentary
  // there for the marker encoding and why the markers become scalar indices in
  // the steady state.
  template<uint32_t L>
  Result run(void)
  {
    using C = Config<L>;
    static_assert(L >= 2 && L <= 127, "outside the signed-byte lane range");

    constexpr uint32_t kTopParity = ((L - 1) & 1u) ^ 1u;

    const __m512i zeros = _mm512_setzero_si512();
    const __m512i ones  = _mm512_set1_epi8(1);
    const __m512i twos  = _mm512_set1_epi8(2);
    const __m512i full  = _mm512_set1_epi8(static_cast<char>(L));

    __m512i w0 = zeros, w1 = zeros;    // segment values, by window
    __m512i m0 = zeros, m1 = zeros;    // one-hot markers, permuted alongside

    uint64_t sweeps = 0;
    int eps = 1;
    bool done = false;

    // --- warm-up: segments still stream in from the initial state ---
    uint32_t parity = 0;
    for (uint32_t step = 0; step < C::kSatCount && !done; step++)
    {
      __m512i curr = parity ? w1 : w0;
      __m512i mcur = parity ? m1 : m0;
      const __mmask64 at_top   = _mm512_cmpeq_epi8_mask(mcur, twos);
      const __mmask64 at_first = _mm512_cmpeq_epi8_mask(mcur, ones);

      if (at_first & _mm512_cmpeq_epi8_mask(curr, full))
      {
        done = true;
        break;
      }
      if (at_first)
      {
        sweeps++;
      }
      if (~_mm512_cmpeq_epi8_mask(curr, zeros) & at_top)
      {
        eps ^= 1;
      }

      const uint32_t p = step;
      const uint8_t value  = (p < L) ? 1 : 0;
      const uint8_t marker = (p == L - 1) ? 2 : ((p == 0) ? 1 : 0);
      __m512i next  = shift_in<C::kLastElem>(parity ? w0 : w1, value);
      __m512i mnext = shift_in<C::kLastElem>(parity ? m0 : m1, marker);

      const __mmask64 eq1 = _mm512_cmpeq_epi8_mask(curr, ones);
      const __mmask64 gt1 = _mm512_cmpgt_epi8_mask(curr, ones);
      __m512i delta = _mm512_add_epi8(_mm512_maskz_mov_epi8(eq1, next),
                                      _mm512_maskz_set1_epi8(gt1, static_cast<char>(-1)));
      delta = _mm512_maskz_mov_epi8(~at_top, delta);

      curr = _mm512_add_epi8(curr, delta);
      next = _mm512_sub_epi8(next, delta);
      if (parity) { w1 = curr; w0 = next; m0 = mnext; }
      else        { w0 = curr; w1 = next; m1 = mnext; }
      parity ^= 1u;
    }

    // --- steady state ---
    constexpr uint32_t kRotPeriod = C::kSatCount / 2;
    uint32_t top_lane = static_cast<uint32_t>(__builtin_ctzll(
      _mm512_cmpeq_epi8_mask((kTopParity == 0) ? m0 : m1, twos)));
    uint32_t first_lane = static_cast<uint32_t>(__builtin_ctzll(
      _mm512_cmpeq_epi8_mask(m1, ones)));

    while (!done)
    {
      // ---- parity 0: w0 is current, w1 is rotated ----
      {
        __m512i curr = w0;
        if constexpr (kTopParity == 0)
        {
          eps ^= static_cast<int>((~_mm512_cmpeq_epi8_mask(curr, zeros) >> top_lane) & 1ULL);
        }
        __m512i next = rot<C::kLastElem>(w1);
        const __mmask64 eq1 = _mm512_cmpeq_epi8_mask(curr, ones);
        const __mmask64 gt1 = _mm512_cmpgt_epi8_mask(curr, ones);
        __m512i delta = _mm512_add_epi8(_mm512_maskz_mov_epi8(eq1, next),
                                        _mm512_maskz_set1_epi8(gt1, static_cast<char>(-1)));
        if constexpr (kTopParity == 0)
        {
          delta = _mm512_maskz_mov_epi8(~(1ULL << top_lane), delta);
        }
        w0 = _mm512_add_epi8(curr, delta);
        w1 = _mm512_sub_epi8(next, delta);
        first_lane = first_lane ? first_lane - 1 : kRotPeriod - 1;
        if constexpr (kTopParity == 1)
        {
          top_lane = top_lane ? top_lane - 1 : kRotPeriod - 1;
        }
      }

      // ---- parity 1: w1 is current, w0 is rotated ----
      {
        __m512i curr = w1;
        if ((_mm512_cmpeq_epi8_mask(curr, full) >> first_lane) & 1ULL)
        {
          done = true;
          break;
        }
        sweeps++;
        if constexpr (kTopParity == 1)
        {
          eps ^= static_cast<int>((~_mm512_cmpeq_epi8_mask(curr, zeros) >> top_lane) & 1ULL);
        }
        __m512i next = rot<C::kLastElem>(w0);
        const __mmask64 eq1 = _mm512_cmpeq_epi8_mask(curr, ones);
        const __mmask64 gt1 = _mm512_cmpgt_epi8_mask(curr, ones);
        __m512i delta = _mm512_add_epi8(_mm512_maskz_mov_epi8(eq1, next),
                                        _mm512_maskz_set1_epi8(gt1, static_cast<char>(-1)));
        if constexpr (kTopParity == 1)
        {
          delta = _mm512_maskz_mov_epi8(~(1ULL << top_lane), delta);
        }
        w1 = _mm512_add_epi8(curr, delta);
        w0 = _mm512_sub_epi8(next, delta);
        if constexpr (kTopParity == 0)
        {
          top_lane = top_lane ? top_lane - 1 : kRotPeriod - 1;
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
} // namespace wide
} // namespace blockless
} // namespace snaperz

#else
#define SNAPERZ_BLOCKLESS_HAVE_AVX512 0
#endif
