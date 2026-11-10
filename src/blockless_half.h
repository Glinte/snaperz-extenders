#pragma once

#include "blockless.h"
#include <limits>
#include <stdexcept>

namespace snaperz { namespace blockless { namespace half {
// Certified 16-coordinate reduction F = P^31 A^435 E, for P=12, L>16.
// C(a,z)=(a,A^31 D(z)) reverses F and fixes its canonical root.
// D places particle j at its position after j-1 inverse tail sweeps.
using State = std::vector<uint32_t>;

inline void sweep(State& a, size_t start = 0)
{
  for (size_t i = start; i + 1 < a.size(); ++i)
  {
    if (a[i] == 1) { a[i] += a[i+1]; a[i+1] = 0; }
    else if (a[i] > 1) { --a[i]; ++a[i+1]; }
  }
}

inline void inverse(State& a, size_t start = 0)
{
  for (size_t j = a.size(); j > start + 1; --j)
  {
    const size_t i = j - 2;
    if (!a[i]) continue;
    if (!a[i+1]) { a[i+1] = a[i] - 1; a[i] = 1; }
    else { ++a[i]; --a[i+1]; }
  }
}

struct Count
{
  Result result;
  uint64_t steps;
  uint64_t rewind;
  bool edge;
};

inline Count run(uint32_t length)
{
  if (length <= 16) throw std::invalid_argument("half-cycle requires L > 16");
  const uint32_t m = length - 16;
  State x(m), y(m), reflected(m), scratch(m);
  x[0] = m;
  uint64_t elapsed = 0, steps = 0, period = 0;
  bool edge;
  for (;;)
  {
    // Leave enough room for endpoint corrections and the physical readout.
    if (elapsed > std::numeric_limits<uint64_t>::max()/4 - m - 466)
      throw std::overflow_error("half-cycle pulse count exceeds uint64_t");
    elapsed += x[0] + uint64_t(466);
    ++steps;
    y = x;
    const uint32_t head = y[0], explicit_steps = std::min(head, m - 1);
    y[0] = 0;
    for (uint32_t i = 0; i < explicit_steps; ++i) { ++y[0]; sweep(y); }
    y.back() += head - explicit_steps;
    for (unsigned i = 0; i < 435; ++i) sweep(y);
    for (unsigned i = 0; i < 31; ++i) sweep(y, 1);

    scratch = y;
    std::fill(reflected.begin(), reflected.end(), 0);
    reflected[0] = y[0];
    const uint32_t mass = m - y[0];
    for (uint32_t rank = 1; rank <= mass; ++rank)
    {
      uint32_t total = 0;
      for (size_t j = 1; j < m; ++j)
      {
        total += scratch[j];
        if (total >= rank) { ++reflected[j]; break; }
      }
      if (rank < mass) inverse(scratch, 1);
    }
    for (unsigned i = 0; i < 31; ++i) sweep(reflected, 1);
    // Edge first: necessary also for a singleton family cycle.
    if (reflected == x) { edge = true; period = 2*elapsed - (m + uint64_t(466)); break; }
    if (reflected == y) { edge = false; period = 2*elapsed + y[0] - m; break; }
    x.swap(y);
  }
  int eps = 0;
  uint64_t rewind = 0;
  if (!(period & 1))
  {
    rewind = elapsed - period/2;
    for (unsigned i = 0; i < 15; ++i) inverse(y, 1);
    State decoded(length);
    decoded[0] = y[0] + 16;
    std::copy(y.begin()+1, y.end(), decoded.begin()+17);
    for (uint64_t i = 0; i < rewind; ++i) inverse(decoded);
    eps = decoded.back() != 0;
  }
  return {{period - (length - 1), eps, period,
           (1 + uint64_t(eps))*period - (length - 1)}, steps, rewind, edge};
}
}}} // namespace snaperz::blockless::half
