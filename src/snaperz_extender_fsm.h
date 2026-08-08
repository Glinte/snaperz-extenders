#pragma once

#include <emmintrin.h> // SSE2
#include <cstring>
#include <cassert>

#include "constants.h"

// An SSE2 implementation for machines without AVX2.
//
// The AVX2 implementation parallelizes *across* pulses: it keeps a window of
// concurrent pulses in flight and advances all of them by one segment per step.
// That is the right thing to do when you have 256-bit integer lanes. This one
// parallelizes *within* a single pulse instead, which needs no AVX2 at all --
// SSE2 is guaranteed on every x86-64 chip -- and is a large win over the scalar
// linked-list fallback.
//
// The trick is a change of coordinates. Instead of storing the segment lengths
// len[k], store their prefix sums
//
//   cum[k] = len[0] + len[1] + ... + len[k],    cum[-1] = 0.
//
// A pulse walks segments left to right, and each segment either does nothing,
// pulls the next segment into itself, or pushes kPushLimit blocks forward. In
// length coordinates every one of those is a *two* cell update (blocks leave
// one segment and arrive in another), which is why the scalar version has to
// chase a linked list. In prefix-sum coordinates all three collapse to a write
// of the single cell cum[k], because moving blocks across the boundary between
// k and k+1 changes exactly that one partial sum:
//
//   len[k] == 0   ->  cum[k] unchanged
//   len[k] == 1   ->  cum[k] = cum[k + 1]     (pull: absorb the next segment)
//   len[k] >  1   ->  cum[k] = cum[k] - 1     (push one block forward)
//
// where len[k] is read as cum[k] minus the *already updated* cum[k - 1].
//
// Note the last line assumes kPushLimit == 1, i.e. period 12. See the guard in
// snaperz_extender.h; other periods use the scalar fallback.
//
// That still looks serial, since cell k needs the updated cell k - 1. But now
// look at what actually crosses the boundary: the carry
//
//   delta[k] = cum'[k] - cum[k]
//
// takes only three values -- 0 (Z, unchanged), -1 (D, pushed), or "jumped to
// cum[k + 1]" (W, pulled) -- and which one fires next depends on the segment
// only through min(len[k], 2), never on how long the segment actually is:
//
//              min(len,0)=0   =1   >=2
//      from Z:            Z    W     D
//      from D:            W    D     D
//      from W:            Z    Z     Z
//
// So a pulse is a three state finite automaton run over the segments, and it
// factors into three passes:
//
//   1. compute the symbols min(len[k], 2)          -- pointwise, vectorized
//   2. run the automaton to get every carry state  -- serial, but table driven
//      at four segments per L1 lookup
//   3. write cum'[k] from the carry state          -- pointwise, vectorized
//
// Pass 3 is pointwise *only* in prefix-sum coordinates; in length coordinates
// it is a two cell update and will not vectorize. That change of coordinates is
// the whole reason this works.
//
// The final segment is left out of the automaton and handled scalar-ly at the
// end, because it is the one place where kLastPushLimit applies rather than
// kPushLimit. It costs at most kLastPushLimit + 1 iterations per pulse.

namespace snaperz
{
  static_assert(std::numeric_limits<len_t>::max() <= std::numeric_limits<uint8_t>::max(),
                "The SSE2 implementation only handles byte-sized segments");
  static_assert(kPushLimit == 1, "The SSE2 implementation assumes a push limit of one");

  // Number of segments, and the total number of blocks distributed among them
  // (the kLength pistons plus the single extended block).
  static constexpr uint32_t kSegCount = kLength + 1;
  static constexpr uint8_t kTotal = static_cast<uint8_t>(kLength + 1);
  // Slack in front of and behind the segments, so that pass 1 may read cum[-1]
  // and the vector passes may over-read and over-write past the last segment.
  static constexpr uint32_t kFront = 16;
  static constexpr uint32_t kBack = 64;
  static constexpr uint32_t kMaskCount = (kSegCount + 2 * 16) / 16;

  struct Extender
  {
    // Prefix sums of the segment lengths. Owns storage starting kFront bytes
    // earlier, which is zeroed so that cum[-1] reads as 0.
    uint8_t* cum;
    // Index of the last non-empty segment, equivalently the first index at
    // which the prefix sum reaches kTotal.
    uint32_t last;
  };

  namespace fsm
  {
    // Carry states. The numeric values are baked into the transition table and
    // into the blends in pass 3.
    static constexpr uint8_t kZ = 0; // cum[k] unchanged
    static constexpr uint8_t kD = 1; // cum[k] -= 1        (pushed a block)
    static constexpr uint8_t kW = 2; // cum[k] = cum[k+1]  (pulled next segment)

    struct Table
    {
      // [carry in][four symbols, packed two bits apiece] -> four carry states,
      // one per byte, the top byte doubling as the carry out.
      uint32_t entries[3][256];
    };

    // Four symbols per lookup keeps the table at 3 KiB, which stays in L1 on
    // any machine. Widening it to eight symbols halves the number of dependent
    // lookups and measured 11.2 s against 15.7 s at length 49 on a Skylake
    // Xeon, but it needs a 1.5 MiB table. That trade only pays off on a chip
    // with a large L2, and this implementation exists precisely for chips old
    // enough to lack AVX2, so the small table is the better default.

    // Flattened, and with the index visibly clamped, so that the compiler can
    // see it is in range while evaluating the table at compile time.
    constexpr uint8_t next_state(uint32_t state, uint32_t symbol)
    {
      const uint8_t transitions[9] = {
        /* from kZ */ kZ, kW, kD,
        /* from kD */ kW, kD, kD,
        /* from kW */ kZ, kZ, kZ,
      };
      const uint32_t index = state * 3 + symbol;
      return transitions[index < 9 ? index : 0];
    }

    constexpr Table build_table()
    {
      Table table{};
      for (uint32_t carry_in = 0; carry_in < 3; carry_in++)
      {
        for (uint32_t index = 0; index < 256; index++)
        {
          uint32_t state = carry_in;
          uint32_t packed = 0;
          for (uint32_t k = 0; k < 4; k++)
          {
            // Symbol k has its low bit in the low nibble of the index and its
            // high bit in the high nibble, which is how the two movemasks in
            // pass 1 lay them out.
            const uint32_t symbol = ((index >> k) & 1) | (((index >> (4 + k)) & 1) << 1);
            state = next_state(state, symbol);
            packed |= state << (8 * k);
          }
          table.entries[carry_in][index] = packed;
        }
      }
      return table;
    }

    static constexpr Table kTable = build_table();

    // Scratch shared by every extender; only live for the duration of a pulse.
    static uint8_t _states[kSegCount + kBack];
    static uint32_t _low_bits[kMaskCount];
    static uint32_t _high_bits[kMaskCount];

    // Runs the automaton over segments [0, count), writing cum in place. The
    // final segment is deliberately excluded; see _simulate_tail.
    inline void _simulate_body(uint8_t* cum, uint32_t count)
    {
      const __m128i _twos = _mm_set1_epi8(2);
      const __m128i _ones = _mm_set1_epi8(1);
      const __m128i _kD = _mm_set1_epi8(kD);
      const __m128i _kW = _mm_set1_epi8(kW);

      // Round up to whole four-segment automaton blocks, and again to whole
      // vectors. The extra segments are harmless: they are cleared out of the
      // state array before pass 3 so that they cannot write anything.
      const uint32_t blocks = (count + 3) / 4;
      const uint32_t vectors = (blocks * 4 + 15) / 16;

      // Pass 1: symbol[k] = min(len[k], 2) = min(cum[k] - cum[k - 1], 2),
      // reduced immediately to a pair of bitmasks holding its two bits.
      for (uint32_t v = 0; v < vectors; v++)
      {
        const __m128i _curr = _mm_loadu_si128((const __m128i*)(cum + v * 16));
        const __m128i _prev = _mm_loadu_si128((const __m128i*)(cum + v * 16 - 1));
        const __m128i _symbol = _mm_min_epu8(_mm_sub_epi8(_curr, _prev), _twos);
        // Shifting a 16-bit lane left by 7 moves bit 0 of each byte into that
        // byte's sign bit, which is what movemask collects; by 6, bit 1.
        _low_bits[v] = (uint32_t)_mm_movemask_epi8(_mm_slli_epi16(_symbol, 7));
        _high_bits[v] = (uint32_t)_mm_movemask_epi8(_mm_slli_epi16(_symbol, 6));
      }

      // Pass 2: the serial part. One L1-resident table lookup advances the
      // automaton by four segments and emits their four carry states.
      uint32_t carry = kZ;
      for (uint32_t b = 0; b < blocks; b++)
      {
        const uint32_t v = b / 4;
        const uint32_t shift = (b % 4) * 4;
        const uint32_t index = ((_low_bits[v] >> shift) & 0xF) |
                               (((_high_bits[v] >> shift) & 0xF) << 4);
        const uint32_t packed = kTable.entries[carry][index];
        std::memcpy(_states + b * 4, &packed, sizeof(packed));
        carry = packed >> 24;
      }
      // Silence the states of the padding segments so pass 3 leaves them alone.
      _mm_storeu_si128((__m128i*)(_states + count), _mm_setzero_si128());

      // Pass 3: cum'[k] follows from the carry state and the old prefix sums.
      for (uint32_t v = 0; v < vectors; v++)
      {
        const __m128i _curr = _mm_loadu_si128((const __m128i*)(cum + v * 16));
        const __m128i _next = _mm_loadu_si128((const __m128i*)(cum + v * 16 + 1));
        const __m128i _state = _mm_loadu_si128((const __m128i*)(_states + v * 16));
        const __m128i _pushed = _mm_cmpeq_epi8(_state, _kD);
        const __m128i _pulled = _mm_cmpeq_epi8(_state, _kW);
        const __m128i _moved = _mm_or_si128(_pushed, _pulled);
        // kD ? curr - 1 : kW ? next : curr
        __m128i _result = _mm_and_si128(_pushed, _mm_sub_epi8(_curr, _ones));
        _result = _mm_or_si128(_result, _mm_and_si128(_pulled, _next));
        _result = _mm_or_si128(_result, _mm_andnot_si128(_moved, _curr));
        _mm_storeu_si128((__m128i*)(cum + v * 16), _result);
      }
    }

    // Handles the segments from the last non-empty one onwards, and returns the
    // index of the last non-empty segment afterwards. This is split out because
    // the last segment holds the extended block, which does not count towards
    // the virtual push limit, so it pushes under kLastPushLimit instead.
    inline uint32_t _simulate_tail(uint8_t* cum, uint32_t last)
    {
      uint32_t k = last;
      uint32_t length = cum[k] - cum[k - 1];
      if (length == 0)
      {
        // The body pulled this segment into its predecessor, so the extender
        // now ends one segment earlier.
        return k - 1;
      }
      while (length > 1)
      {
        const uint32_t pushed = std::min<uint32_t>(kLastPushLimit, length - 1);
        cum[k] -= pushed;
        // Everything past k still holds the full block count, so the segment
        // that just received the blocks has exactly that many.
        k++;
        length = pushed;
      }
      return k;
    }
  } // namespace fsm

  Extender create()
  {
    Extender extender;
    uint8_t* storage = new uint8_t[kFront + kSegCount + kBack];
    std::memset(storage, 0, kFront);
    extender.cum = storage + kFront;
    // Fully extended: every segment holds one block, so the prefix sums count
    // up from one. Past the end they stay saturated at the total block count,
    // which is what makes the padding segments read as empty.
    for (uint32_t i = 0; i < kSegCount; i++)
    {
      extender.cum[i] = static_cast<uint8_t>(i + 1);
    }
    std::memset(extender.cum + kSegCount, kTotal, kBack);
    extender.last = kLength;
    return extender;
  }

  void destroy(Extender& extender)
  {
    delete[] (extender.cum - kFront);
    extender.cum = nullptr;
  }

  void simulate_pulse(Extender& extender)
  {
    assert(extender.last > 0 && "the extender has already finished");
    fsm::_simulate_body(extender.cum, extender.last);
    extender.last = fsm::_simulate_tail(extender.cum, extender.last);
  }

  bool equals(const Extender& lhs, const Extender& rhs)
  {
    return std::memcmp(lhs.cum, rhs.cum, kSegCount) == 0;
  }

  bool finished(const Extender& extender)
  {
    // Every block sits in the first segment.
    return extender.cum[0] == kTotal;
  }
}
