"""Round-28: the arithmetic hiding in the sweep's gate factorization.

The round-22 result ord(g_i) = lcm(1, ..., n-i) is the one statement in the
ledger whose *shape* does not match the rest of it: an arithmetic object
(lcm(1..k) = e^{psi(k)}, Chebyshev) sitting inside a combinatorial dynamical
system.  A_n factors as g_0 g_1 ... g_{n-2}, so the sweep is a product of n-1
permutations with completely known, arithmetically clean orders.

This measures what that factorization actually buys:

  1. the exact cycle-length multiset of each gate g_i (is it {1, ..., n-i}?);
  2. ord(A_n), the lcm of all A_n cycle lengths, and its factorization --
     do the large primes in p_L come from the gates or only from the product?
  3. whether the primes dividing p_L are bounded by anything gate-related.

usage: python3 r28_gatearith.py [nmax]
"""
import sys
from collections import Counter
from math import lcm

import core


def gate(x, i):
    x = list(x)
    a = x[i]
    if a == 1:
        x[i] = 1 + x[i + 1]
        x[i + 1] = 0
    elif a >= 2:
        x[i] = a - 1
        x[i + 1] += 1
    return tuple(x)


def cycle_lengths(states, step):
    """Multiset of cycle lengths of a permutation given as a step function."""
    seen, out = set(), Counter()
    for x in states:
        if x in seen:
            continue
        c, y = 0, x
        while y not in seen:
            seen.add(y)
            y = step(y)
            c += 1
        out[c] += 1
    return out


def factor(m):
    f, d = Counter(), 2
    while d * d <= m:
        while m % d == 0:
            f[d] += 1
            m //= d
        d += 1
    if m > 1:
        f[m] += 1
    return f


def fmt(f):
    return " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items()))


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    print("gate cycle structure: is the length set exactly {1, ..., n-i}?")
    for n in range(2, min(nmax, 10) + 1):
        states = core.catalan_states(n)
        row = []
        for i in range(n - 1):
            cl = cycle_lengths(states, lambda x, i=i: gate(x, i))
            want = set(range(1, n - i + 1))
            got = set(cl)
            row.append("ok" if got == want else f"i={i}:{sorted(got)}")
        print(f"  n={n:2d}: " + " ".join(row))

    print("\norder of the whole sweep, and where its primes come from:")
    print(f"  {'n':>3} {'|C_n|':>8} {'p_n':>8} {'cycles':>7} {'max':>7} "
          f"{'ord(A_n)':>22}  factorization of ord(A_n)")
    for n in range(2, nmax + 1):
        states = core.catalan_states(n)
        cl = cycle_lengths(states, core.A)
        order = 1
        for c in cl:
            order = lcm(order, c)
        R = (n,) + (0,) * (n - 1)
        p, y = 0, R
        while True:
            y = core.A(y)
            p += 1
            if y == R:
                break
        mx = max(cl)
        os_ = str(order) if order < 10 ** 20 else f"~10^{len(str(order)) - 1}"
        print(f"  {n:3d} {len(states):8d} {p:8d} {sum(cl.values()):7d} {mx:7d} "
              f"{os_:>22}  {fmt(factor(order)) if order < 10**18 else '(large)'}")

    print("\nprimes of p_n vs n, and vs lcm(1..n):")
    for n in range(4, nmax + 1):
        R = (n,) + (0,) * (n - 1)
        p, y = 0, R
        while True:
            y = core.A(y)
            p += 1
            if y == R:
                break
        f = factor(p)
        big = max(f) if f else 1
        L = lcm(*range(1, n + 1))
        print(f"  n={n:2d} p_n={p:8d} = {fmt(f):28s} largest prime {big:6d} "
              f"{'<= n' if big <= n else '> n'}   p_n | lcm(1..n)? "
              f"{'yes' if L % p == 0 else 'no'}")


if __name__ == "__main__":
    main()
