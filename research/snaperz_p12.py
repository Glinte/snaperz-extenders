#!/usr/bin/env python3
"""Exact pulse counts for P=12 Snaperz extenders.

This collects every algorithmic result from the derivation:

1.  A state is the tuple of air-separated segment lengths ``(a_0, ..., a_{n-1})``
    with ``sum(a) == n``.  The reachable states are exactly the *Catalan states*
    ``a_0 + ... + a_k >= k + 1``, i.e. Dyck paths of semilength ``n``.
2.  Deleting the extender block leaves the *ordinary blockless sweep* ``A``,
    which is a permutation of the Catalan states.
3.  ``E_L = (1, ..., 1)`` (extended) and ``R_L = (L, 0, ..., 0)`` (retracted)
    lie on one ``A``-cycle, of length ``p_L``.  Since ``A^(L-1)(R_L) = E_L``,

        B_L = p_L - (L - 1)

    is the first-hitting time of ``R_L`` from ``E_L``.
4.  Re-inserting the block adds one fibre coordinate ``r``.  The good-pair
    invariant tracks the real block together with an auxiliary block; exactly
    one of them retracts on the first trip.  The winner bit is

        eps_L = parity of the number of role swaps around the cycle,

    and a role swap happens exactly on a sweep whose final-tail height is 1.
5.  T_L = (1 + eps_L) * p_L - (L - 1).
6.  Checkpoints are the states with ``a_0 == 1``.  The induced *checkpoint map*
    ``H`` steps between them, giving

        q_m = c_{m+1}(1)          (number of checkpoints in one period)
        S_m = p_{m+1} - q_m
        B_L = q_{L-1} + S_{L-1} - (L - 1).
7.  Nested first-return maps ``F_k`` (the section of ``H`` below a stem of ``k``
    leading ones) renormalise the dynamics.  ``F_k`` has a *constant* H-return
    time for k in {1, 3, 7, 8, 16}, of 1, 2, 8, 8 and 112 H-steps.  The last one
    is the ``H^112(1^15 u) = 1^15 J(u)`` law, and it forces

        q_m = 112 * j_{m-15}          for m >= 15.

Performance note: the nested ``F_k`` recursion telescopes exactly to the plain
``A``-orbit of ``E_L`` -- memoising it never produces a single cache hit, since
each induced state is visited once.  So the hot path here is a single flat
sweep loop, optionally compiled to C.  The 112x renormalisation speedup is *not*
available: it would require evaluating ``J`` faster than the 112 H-steps that
define it, and ``J`` has no known closed form (it is not a power of ``A``).
``renormalised_period`` below still reports ``j``, and ``verify_sections``
checks the constant-return-time laws.

Usage:
    python3 snaperz_p12.py 30 31
    python3 snaperz_p12.py --table 40
    python3 snaperz_p12.py --self-test
"""

import argparse
import ctypes
import hashlib
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass

type State = tuple[int, ...]

# H-return times of the induced section maps F_k that are independent of the
# suffix.  F_1 = H, F_16 = J.  Verified exhaustively for all Catalan suffixes
# up to size 6 by verify_sections().
CONSTANT_TIME_SECTIONS: dict[int, int] = {1: 1, 3: 2, 7: 8, 8: 8, 16: 112}

RENORMALISATION_STEM = 16
RENORMALISATION_TIME = CONSTANT_TIME_SECTIONS[RENORMALISATION_STEM]


# --------------------------------------------------------------------------- #
# 1. State space
# --------------------------------------------------------------------------- #


def is_catalan_state(state: State) -> bool:
    """Return whether ``state`` satisfies the prefix condition and the total."""
    running = 0
    for index, value in enumerate(state):
        running += value
        if running < index + 1:
            return False
    return running == len(state)


def to_dyck_word(state: State) -> str:
    """Render a Catalan state as balanced parentheses, ``(1, 1, 1) -> ()()()``."""
    return "".join("(" * value + ")" for value in state)


def extended_state(length: int) -> State:
    """Return ``E_n = (1, ..., 1)``."""
    return (1,) * length


def retracted_state(length: int) -> State:
    """Return ``R_n = (n, 0, ..., 0)``."""
    return (length,) + (0,) * (length - 1)


# --------------------------------------------------------------------------- #
# 2. The ordinary blockless sweep A
# --------------------------------------------------------------------------- #


def ordinary_sweep(state: State) -> State:
    """Apply one ordinary blockless P=12 sweep ``A``.

    Adjacent pairs are processed left to right, and later pairs see the updates
    made by earlier ones::

        (0, y) -> (0, y)
        (1, y) -> (1 + y, 0)
        (x, y) -> (x - 1, y + 1)    for x > 1

    In Dyck-word language each closing bracket either steps one place left, or
    (when only one opening bracket precedes it) hops right to just before the
    next closing bracket.  This is the plane-tree whirling action.
    """
    result = list(state)
    for index in range(len(result) - 1):
        current = result[index]
        if current == 0:
            continue
        following = result[index + 1]
        if current == 1:
            result[index] = 1 + following
            result[index + 1] = 0
        else:
            result[index] = current - 1
            result[index + 1] = following + 1
    return tuple(result)


def final_tail_height(state: State) -> int:
    """Return the final-descent height ``h``: ``n - j`` for the last nonempty ``j``."""
    for index in range(len(state) - 1, -1, -1):
        if state[index] > 0:
            return len(state) - index
    raise ValueError("A piston state cannot be entirely empty.")


def checkpoint_step(state: State) -> tuple[State, int]:
    """Apply the checkpoint map ``H`` to the tail of a checkpoint ``(1, state)``.

    Returns the new tail and the number of ``A``-sweeps consumed, which is
    ``state[0] + 1``: the leading segment jumps to ``1 + state[0]`` and then
    counts back down to 1.
    """
    sweeps = state[0] + 1
    full = (1,) + state
    for _ in range(sweeps):
        full = ordinary_sweep(full)
    if full[0] != 1:
        raise AssertionError("The checkpoint did not return as expected.")
    return full[1:], sweeps


# --------------------------------------------------------------------------- #
# 3. The orbit scan (hot loop)
# --------------------------------------------------------------------------- #

_C_SOURCE = r"""
#include <stdint.h>
#include <stdlib.h>

int orbit_data(int L, uint64_t *out_q, uint64_t *out_p, int *out_parity) {
    uint8_t *s = (uint8_t *)malloc((size_t)L + 1);
    if (!s) return -1;
    for (int i = 0; i < L; i++) s[i] = 1;

    const int top = L - 1;
    int last = top, parity = 0;
    uint64_t b = 0, q = 0;

    for (;;) {
        if (last == top) parity ^= 1;   /* final-tail height 1: roles swap */
        if (s[0] == 1) q++;             /* checkpoint */

        int stop = last + 1;
        if (stop > top) stop = top;
        for (int i = 0; i < stop; i++) {
            int c = s[i];
            if (!c) continue;
            int n = s[i + 1];
            if (c == 1) { s[i] = (uint8_t)(1 + n); s[i + 1] = 0; }
            else        { s[i] = (uint8_t)(c - 1); s[i + 1] = (uint8_t)(n + 1); }
        }
        b++;

        if (last < top && s[last + 1]) last++;
        while (!s[last]) last--;
        if (s[0] == L) break;           /* reached R_L */
    }

    free(s);
    *out_q = q; *out_p = b + (uint64_t)(L - 1); *out_parity = parity;
    return 0;
}
"""

_c_orbit_data: "ctypes.CDLL | None | bool" = False


def _load_c_core() -> object | None:
    """Compile and load the C sweep loop, or return None if that is impossible.

    The build is cached in the temp directory under a hash of the source, so it
    happens at most once per machine per version of this file.
    """
    global _c_orbit_data
    if _c_orbit_data is not False:
        return _c_orbit_data

    _c_orbit_data = None
    digest = hashlib.sha256(_C_SOURCE.encode()).hexdigest()[:16]
    library_path = os.path.join(tempfile.gettempdir(), f"snaperz_p12_{digest}.so")
    try:
        if not os.path.exists(library_path):
            source_path = library_path + ".c"
            with open(source_path, "w") as handle:
                handle.write(_C_SOURCE)
            subprocess.run(
                ["cc", "-O3", "-shared", "-fPIC", "-o", library_path, source_path],
                check=True,
                capture_output=True,
            )
        library = ctypes.CDLL(library_path)
        library.orbit_data.argtypes = [
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_uint64),
            ctypes.POINTER(ctypes.c_uint64),
            ctypes.POINTER(ctypes.c_int),
        ]
        library.orbit_data.restype = ctypes.c_int
        _c_orbit_data = library
    except (OSError, subprocess.SubprocessError):
        _c_orbit_data = None
    return _c_orbit_data


def _orbit_data_python(length: int) -> tuple[int, int, int]:
    """Pure-Python fallback for :func:`orbit_data`."""
    state = [1] * length
    top = length - 1
    last = top
    hits = 0
    checkpoints = 0
    parity = 0

    while True:
        if last == top:
            parity ^= 1
        if state[0] == 1:
            checkpoints += 1

        # Everything past the last nonempty segment is fixed by the sweep, and
        # the one segment the sweep can spill into is a no-op on the next pair.
        stop = last + 1
        if stop > top:
            stop = top
        for index in range(stop):
            current = state[index]
            if current:
                following = state[index + 1]
                if current == 1:
                    state[index] = 1 + following
                    state[index + 1] = 0
                else:
                    state[index] = current - 1
                    state[index + 1] = following + 1
        hits += 1

        if last < top and state[last + 1]:
            last += 1
        while state[last] == 0:
            last -= 1

        # Stopping at R_L costs one integer compare instead of an O(L) scan for
        # E_L: no other state of sum L over L segments can lead with L.
        if state[0] == length:
            return checkpoints, hits + top, parity


def orbit_data(length: int, *, use_c: bool = True) -> tuple[int, int, int]:
    """Scan one full ``A``-orbit of ``E_L`` and return ``(q, p_L, eps_L)``.

    ``q`` counts the checkpoints (states with a leading segment of 1) and is the
    checkpoint period ``q_{L-1}``; ``p_L`` is the ordinary period; ``eps_L`` is
    the good-pair winner bit.

    The winner bit is accumulated as the parity of the sweeps whose *starting*
    state has final-tail height 1, rather than whose *ending* state does.  A
    role swap actually happens on the latter, but around a closed cycle every
    state is the start of one sweep and the end of another, so the two parities
    agree -- and the starting test is just ``last == top``.

    The scan stops at ``R_L`` after ``B_L`` sweeps rather than closing the loop
    at ``E_L``.  The remaining ``L - 1`` states are ``A^j(R_L) = (L-j, 1^j, 0^..)``
    explicitly, and they contribute no checkpoint (their leading segments run
    ``L`` down to ``2``) and no parity flip (their final segment is empty), so
    ``q`` and ``eps`` are already complete and ``p_L = B_L + (L - 1)``.
    """
    if length < 1:
        raise ValueError("The extender length must be at least 1.")
    if length == 1:
        return 1, 1, 0

    if use_c:
        library = _load_c_core()
        if library is not None:
            checkpoints = ctypes.c_uint64()
            sweeps = ctypes.c_uint64()
            parity = ctypes.c_int()
            if library.orbit_data(length, checkpoints, sweeps, parity) == 0:
                return checkpoints.value, sweeps.value, parity.value

    return _orbit_data_python(length)


# --------------------------------------------------------------------------- #
# 4. The report
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ExtenderReport:
    """Every derived quantity for one P=12 extender length."""

    extender_length: int          # L
    checkpoint_period: int        # q_{L-1} = c_L(1)
    checkpoint_weight: int        # S_{L-1} = p_L - q_{L-1}
    ordinary_period: int          # p_L
    first_hitting_time: int       # B_L = p_L - (L - 1)
    winner_bit: int               # eps_L
    total_pulses: int             # T_L

    @property
    def winner(self) -> str:
        return "real block" if self.winner_bit == 0 else "auxiliary block"

    @property
    def renormalised_period(self) -> int | None:
        """``j_{L-16}`` from ``q_{L-1} = 112 * j_{L-16}``, or None below the stem."""
        if self.extender_length - 1 < RENORMALISATION_STEM - 1:
            return None
        divisor, remainder = divmod(self.checkpoint_period, RENORMALISATION_TIME)
        return divisor if remainder == 0 else None


def analyse(extender_length: int, *, use_c: bool = True) -> ExtenderReport:
    """Compute the full report for an extender with ``extender_length`` pistons."""
    checkpoints, sweeps, parity = orbit_data(extender_length, use_c=use_c)
    tail = extender_length - 1
    return ExtenderReport(
        extender_length=extender_length,
        checkpoint_period=checkpoints,
        checkpoint_weight=sweeps - checkpoints,
        ordinary_period=sweeps,
        first_hitting_time=sweeps - tail,
        winner_bit=parity,
        total_pulses=(1 + parity) * sweeps - tail,
    )


def print_report(report: ExtenderReport) -> None:
    """Print one report in the layout of the original script."""
    print(f"L               = {report.extender_length}")
    print(f"q_(L-1)         = {report.checkpoint_period:,}")
    print(f"S_(L-1)         = {report.checkpoint_weight:,}")
    print(f"p_L             = {report.ordinary_period:,}")
    print(f"B_L             = {report.first_hitting_time:,}")
    print(f"first winner    = {report.winner}")
    print(f"epsilon_L       = {report.winner_bit}")
    print(f"T_L             = {report.total_pulses:,}")
    renormalised = report.renormalised_period
    if renormalised is not None:
        print(f"j_(L-16)        = {renormalised:,}    (q = {RENORMALISATION_TIME} * j)")


# --------------------------------------------------------------------------- #
# 5. Renormalisation: the induced section maps F_k
# --------------------------------------------------------------------------- #


def catalan_states(size: int) -> list[State]:
    """Enumerate every Catalan state of the given size."""
    states: list[State] = []

    def extend(prefix: list[int], total: int) -> None:
        index = len(prefix)
        if index == size:
            if total == size:
                states.append(tuple(prefix))
            return
        for value in range(size - total + 1):
            if total + value >= index + 1:
                prefix.append(value)
                extend(prefix, total + value)
                prefix.pop()

    extend([], 0)
    return states


def induced_section(stem_length: int, suffix: State) -> tuple[State, int, int]:
    """Evaluate ``F_k(suffix)``, the first return of ``H`` to the stem ``1^k``.

    Returns the image suffix, the number of ``H`` steps, and the number of
    ``A`` sweeps.  ``F_1`` is ``H`` itself and ``F_16`` is ``J``.
    """
    stem = (1,) * stem_length
    state = stem + suffix
    h_steps = 0
    a_sweeps = 0
    while True:
        state = ordinary_sweep(state)
        a_sweeps += 1
        if state[0] == 1:
            h_steps += 1
            if state[:stem_length] == stem:
                return state[stem_length:], h_steps, a_sweeps


def verify_parity_quotient(half_size: int) -> tuple[int, int]:
    """Check the parity-phase invariants on the unary orbit of size ``2n + 1``.

    At every second checkpoint the state is ``(1, x_1, y_1, ..., x_n, y_n)`` with
    each pair either ``(odd, odd)`` or ``(0, even)``, the coarse vector
    ``z_i = (x_i + y_i) / 2`` is itself a Catalan state of size ``n``, and the two
    checkpoints together consume exactly ``2 * z_1 + 4`` ordinary sweeps.

    These do not hold for arbitrary Catalan states -- only on the reachable
    parity phase -- which makes them a sharp assertion on the orbit.  Returns
    ``(checkpoint period, number of coarse steps checked)``.
    """
    size = 2 * half_size + 1
    start = extended_state(size)
    state = start
    orbit: list[State] = []
    weights: list[int] = []
    while True:
        orbit.append(state)
        state, sweeps = checkpoint_step(state)
        weights.append(sweeps)
        if state == start:
            break

    period = len(orbit)
    if period % 2:
        raise AssertionError(f"q_{size} = {period} is odd; the quotient is undefined.")

    for index in range(0, period, 2):
        state = orbit[index]
        if state[0] != 1:
            raise AssertionError(f"Checkpoint {index} does not lead with 1: {state}.")
        pairs = [(state[2 * i + 1], state[2 * i + 2]) for i in range(half_size)]
        for x, y in pairs:
            if not ((x % 2 and y % 2) or (x == 0 and y % 2 == 0)):
                raise AssertionError(f"Pair {(x, y)} at checkpoint {index} breaks the parity form.")
        coarse = tuple((x + y) // 2 for x, y in pairs)
        if not is_catalan_state(coarse):
            raise AssertionError(f"Coarse state {coarse} is not Catalan.")
        if weights[index] + weights[index + 1] != 2 * coarse[0] + 4:
            raise AssertionError(
                f"Weight {weights[index] + weights[index + 1]} != 2*{coarse[0]}+4 at checkpoint {index}."
            )

    return period, period // 2


def verify_sections(max_suffix_size: int = 5) -> dict[int, int | None]:
    """Check which ``F_k`` have a suffix-independent H-return time.

    Returns a mapping from stem length to the constant return time, or None
    where the return time depends on the suffix.
    """
    results: dict[int, int | None] = {}
    for stem_length in sorted(CONSTANT_TIME_SECTIONS):
        times = {
            induced_section(stem_length, suffix)[1]
            for size in range(max_suffix_size + 1)
            for suffix in catalan_states(size)
        }
        results[stem_length] = times.pop() if len(times) == 1 else None
    return results


# --------------------------------------------------------------------------- #
# 6. Ground truth: simulate the real extender
# --------------------------------------------------------------------------- #


def full_pulse(state: State) -> State:
    """Apply one complete P=12 pulse to the real extender, block included."""
    result = list(state)
    size = len(result)
    for index in range(size - 1):
        current = result[index]
        if current == 0:
            continue
        if all(result[j] == 0 for j in range(index + 1, size)):
            # The block segment: the next segment is empty by definition.
            if current == 2:
                result[index] = 1
                result[index + 1] = 1
            elif current > 2:
                result[index] = current - 2
                result[index + 1] = 2
        elif current == 1:
            result[index] = 1 + result[index + 1]
            result[index + 1] = 0
        else:
            result[index] = current - 1
            result[index + 1] += 1
    return tuple(result)


def brute_force_pulse_count(extender_length: int, limit: int = 10**7) -> int:
    """Count pulses to retraction by simulating the real extender directly."""
    size = extender_length + 1
    state = extended_state(size)
    target = retracted_state(size)
    for pulses in range(1, limit + 1):
        state = full_pulse(state)
        if state == target:
            return pulses
    raise RuntimeError(f"No retraction within {limit:,} pulses.")


def self_test(max_brute_force: int = 22) -> None:
    """Check the fast path against brute force and the renormalisation laws."""
    for length in range(1, max_brute_force + 1):
        expected = brute_force_pulse_count(length)
        report = analyse(length)
        assert report.total_pulses == expected, (length, report.total_pulses, expected)
        assert report.ordinary_period == report.first_hitting_time + length - 1
        assert report.checkpoint_period + report.checkpoint_weight == report.ordinary_period
    print(f"T_L matches direct extender simulation for L = 1..{max_brute_force}")

    if analyse(1, use_c=False) != analyse(1) or any(
        analyse(length, use_c=False) != analyse(length) for length in range(2, 25)
    ):
        raise AssertionError("The C core disagrees with the Python fallback.")
    print("C core matches the pure-Python fallback for L = 1..24")

    sections = verify_sections()
    assert sections == CONSTANT_TIME_SECTIONS, sections
    print(f"constant-time induced sections F_k confirmed: {sections}")

    checked = [verify_parity_quotient(n) for n in range(1, 10)]
    steps = sum(coarse for _, coarse in checked)
    print(f"parity-phase form, coarse Catalan state and the 2*z_1+4 weight hold over {steps} coarse steps")


# --------------------------------------------------------------------------- #
# 7. Command line
# --------------------------------------------------------------------------- #


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("lengths", type=int, nargs="*", help="Extender lengths L.")
    parser.add_argument("--table", type=int, metavar="MAX_L", help="Print a table up to MAX_L.")
    parser.add_argument("--self-test", action="store_true", help="Verify against brute force.")
    parser.add_argument("--no-c", action="store_true", help="Force the pure-Python core.")
    args = parser.parse_args()

    use_c = not args.no_c

    if args.self_test:
        self_test()
        return

    if args.table:
        print(f"{'L':>4} {'q_(L-1)':>10} {'S_(L-1)':>12} {'p_L':>12} {'B_L':>12} {'eps':>4} {'T_L':>13}")
        for length in range(1, args.table + 1):
            report = analyse(length, use_c=use_c)
            print(
                f"{report.extender_length:>4} {report.checkpoint_period:>10,} "
                f"{report.checkpoint_weight:>12,} {report.ordinary_period:>12,} "
                f"{report.first_hitting_time:>12,} {report.winner_bit:>4} "
                f"{report.total_pulses:>13,}"
            )
        return

    lengths = args.lengths or [30, 31]
    for position, length in enumerate(lengths):
        if position:
            print()
        print_report(analyse(length, use_c=use_c))


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
