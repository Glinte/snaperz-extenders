#include <iostream>
#include <chrono>
#include <cstdint>
#include <iomanip>

#include "blockless.h"
#include "constants.h"

// Counts the pulses of a P=12 snaperz extender without simulating the extender
// block, using the good-pair reduction:
//
//     T_L = (1 + eps_L) * p_L - (L - 1),   p_L = B_L + (L - 1)
//
// where B_L is the number of blockless sweeps from E_L to R_L and eps_L is the
// parity of the sweeps whose starting state has a non-empty top segment. When
// eps_L is 1 this reaches the same answer in half the sweeps that simulating
// the extender directly takes; when it is 0 the two are the same amount of
// work. It is never more.
//
// This is deliberately a separate program rather than a mode of the simulator.
// The simulator arrives at T_L by simulating the real extender, and this
// arrives at it from the theory; keeping them apart keeps them independent, so
// that they can be used to check each other.
int main()
{
  std::cout
    << "Running blockless " << kLength << " extender, "
    << kPeriod << " tick period." << std::endl;

  const auto start = std::chrono::steady_clock::now();
  const snaperz::blockless::Result result = snaperz::blockless::run();
  const double seconds = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - start).count();

  std::cout
    << "B_L     = " << result.sweeps << '\n'
    << "eps_L   = " << result.eps
    << (result.eps ? "  (auxiliary block wins, so T_L is about 2 * p_L)"
                   : "  (real block wins)") << '\n'
    << "p_L     = " << result.p << '\n'
    << "T_L     = " << result.total << " pulses" << '\n'
    << "Done in " << std::fixed << std::setprecision(2) << seconds << "s"
    << std::endl;
  return 0;
}
