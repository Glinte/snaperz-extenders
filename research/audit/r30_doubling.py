"""Round-30: does the doubling map on the canonical cycle have local structure?

The unique-marker no-go leaves one algorithmic shape open: an injective
representation with cheap random access.  The canonical such object is

    D : x_t -> x_{2t mod p}      (canonical cycle rooted at R_L)

If D were poly(L)-computable from the state alone, double-and-add would give
T_L in poly(L) log p -- the subexponential algorithm.  D is injective, so the
no-go does not forbid it.  Nothing in 29 rounds has measured it.

Diagnostic: prefix determinacy.  For prefix length m, what fraction of cycle
states have output coordinate j determined by input prefix x_0..x_{m-1}
(i.e. all states sharing that prefix agree on the output coordinate)?

Calibrations:
  A      : one sweep, prefix-causal with lag 1 -- determined at m = j+2;
  A^-1   : anti-causal, should look bad in prefix terms;
  REFL   : t -> -t, the reversor's on-cycle action.  Cheap via the ECO
           recursion yet plausibly prefix-nonlocal: calibrates "nonlocal in
           this metric but cheap anyway", so a bad D curve is evidence,
           not proof, of death;
  RAND   : a uniformly shuffled cycle permutation -- the no-structure null.
           If D tracks RAND, D carries no exploitable local structure.

usage: python3 r30_doubling.py [L ...]     (defaults 14 16 18 20 21 22)
"""
import random
import sys

import core


def cycle(L):
    R = (L,) + (0,) * (L - 1)
    states, x = [], R
    while True:
        states.append(x)
        x = core.A(x)
        if x == R:
            break
    return states


def determinacy(states, out, j, mmax):
    """fraction of states whose output coord j is forced by prefix m, per m."""
    p = len(states)
    res = []
    for m in range(1, mmax + 1):
        groups = {}
        for t, x in enumerate(states):
            groups.setdefault(x[:m], set()).add(out[t][j])
        forced = sum(1 for x in states if len(groups[x[:m]]) == 1)
        res.append(forced / p)
    return res


def m_star(curve, frac):
    for m, v in enumerate(curve, 1):
        if v >= frac:
            return m
    return None


def main():
    Ls = [int(a) for a in sys.argv[1:]] or [14, 16, 18, 20, 21, 22]
    random.seed(0)
    print(f"{'L':>3} {'p':>6}  map    m*: coord0  coord1  coord2   (prefix length forcing 99% of states)")
    for L in Ls:
        st = cycle(L)
        p = len(st)
        maps = {
            "A":    [st[(t + 1) % p] for t in range(p)],
            "A^-1": [st[(t - 1) % p] for t in range(p)],
            "REFL": [st[(-t) % p] for t in range(p)],
            "D":    [st[(2 * t) % p] for t in range(p)],
            "T3":   [st[(3 * t) % p] for t in range(p)],
            "RAND": random.sample(st, p),
        }
        for name, out in maps.items():
            ms = [m_star(determinacy(st, out, j, L), 0.99) for j in (0, 1, 2)]
            ms = [str(m) if m else ">L" for m in ms]
            print(f"{L:>3} {p:>6}  {name:<5}  {ms[0]:>9} {ms[1]:>7} {ms[2]:>7}")
        print()


if __name__ == "__main__":
    main()
