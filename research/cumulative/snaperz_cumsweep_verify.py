"""Does the cumulative representation buy anything for the P=12 sweep?"""

import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from snaperz_p12 import ordinary_sweep, catalan_states, extended_state  # noqa: E402


def to_cum(a):
    s, acc = [], 0
    for x in a:
        acc += x
        s.append(acc)
    return s


def from_cum(s):
    prev, a = 0, []
    for x in s:
        a.append(x - prev)
        prev = x
    return tuple(a)


def sweep_cum(a):
    """The sweep written directly on s_i = sum_{j<=i} a_j.  One write per gate."""
    s = to_cum(a)
    prev = 0
    for i in range(len(s) - 1):
        cur = s[i]
        if cur == prev:  # a_i = 0
            pass
        elif cur == prev + 1:  # a_i = 1: wrap
            cur = s[i + 1]
            s[i] = cur
        else:  # a_i >= 2
            cur -= 1
            s[i] = cur
        prev = cur
    return from_cum(s)


# The carry between gates is delta_i = s_i' - s_i, and it only ever takes three
# values: 0 (Z), -1 (D), or "jumped to s_{i+1}" (W).  Which one fires depends on
# the incoming carry and on min(a_i, 2) only -- never on the magnitude of a_i.
Z, D, W = 0, 1, 2
DELTA = {  # (state_in, min(a_i,2)) -> state_out
    (Z, 0): Z, (Z, 1): W, (Z, 2): D,
    (D, 0): W, (D, 1): D, (D, 2): D,
    (W, 0): Z, (W, 1): Z, (W, 2): Z,
}


def sweep_automaton(a):
    """Sweep as: symbols -> 3-state carry scan -> pointwise apply."""
    s = to_cum(a)
    n = len(s) - 1
    sym = [min(s[i] - (s[i - 1] if i else 0), 2) for i in range(n)]  # pass 1
    st, states = Z, []
    for i in range(n):  # pass 2: the only serial part
        st = DELTA[(st, sym[i])]
        states.append(st)
    out = list(s)
    for i in range(n):  # pass 3: pointwise, from OLD s only
        if states[i] == D:
            out[i] = s[i] - 1
        elif states[i] == W:
            out[i] = s[i + 1]
    return from_cum(out)


def main():
    bad = 0
    checked = 0
    for L in range(2, 13):
        for st in catalan_states(L):
            want = ordinary_sweep(st)
            for name, fn in (("cum", sweep_cum), ("auto", sweep_automaton)):
                got = fn(st)
                checked += 1
                if got != want:
                    bad += 1
                    if bad < 6:
                        print(f"MISMATCH {name} L={L} {st} -> {got} want {want}")
    # also along real orbits, where the interesting states actually live
    for L in (17, 31, 41):
        st = extended_state(L)
        for _ in range(20000):
            want = ordinary_sweep(st)
            for name, fn in (("cum", sweep_cum), ("auto", sweep_automaton)):
                got = fn(st)
                checked += 1
                if got != want:
                    bad += 1
                    if bad < 6:
                        print(f"ORBIT MISMATCH {name} L={L} {st} -> {got} want {want}")
            st = want
    print(f"checked={checked} mismatches={bad}")

    # How much of a sweep is actually a no-op (free in cumulative coords)?
    for L in (31, 41, 49):
        st = extended_state(L)
        zeros = tot = 0
        for _ in range(200000):
            zeros += sum(1 for x in st[: L - 1] if x == 0)
            tot += L - 1
            st = ordinary_sweep(st)
        print(f"L={L}: a_i==0 fraction over orbit prefix = {zeros / tot:.3f}")


main()
