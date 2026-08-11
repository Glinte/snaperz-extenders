"""Round-20/23 audit: the even sector as 2-Dyck states, K = A^2, the no-flat
theorem, and the claimed one-bit-carry sequential kernel for the half-sweep.

usage: python3 r23_even.py [nmax]
"""
import sys
from math import comb

import core


def enc_pair(u, v):
    if u == 0:
        assert v % 2 == 0
        return (0, v // 2)
    assert u % 2 == 1 and v % 2 == 1
    return ((u + 1) // 2, (v - 1) // 2)


def enc_aligned(x):
    n = len(x) // 2
    b = []
    for j in range(n):
        b.extend(enc_pair(x[2 * j], x[2 * j + 1]))
    return tuple(b)


def enc_staggered(y):
    """(q, c_0, d_0, ..., c_{n-2}, d_{n-2}, 0)"""
    n = len(y) // 2
    assert y[0] > 0 and y[0] % 2 == 0 and y[2 * n - 1] == 0
    out = [y[0] // 2]
    for j in range(n - 1):
        out.extend(enc_pair(y[2 * j + 1], y[2 * j + 2]))
    out.append(0)
    return tuple(out)


def is_2dyck(b, n):
    if sum(b) != n or len(b) != 2 * n:
        return False
    s = 0
    for k, v in enumerate(b):
        s += v
        if 2 * s < k + 1:
            return False
    return True


def kappa(e, c, d):
    """The claimed one-bit-carry kernel."""
    if e == 0:
        if c == 0:
            return (0, d), (0 if d == 0 else 1)
        if c == 1:
            return (d + 1, 0), 0
        return (c - 1, d + 1), 1
    if c == 0:
        return (d, 0), 0
    return (c - 1, d + 1), 1


def hgt(x):
    m = len(x) - 1
    while m >= 0 and x[m] == 0:
        m -= 1
    return len(x) - m


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for n in range(1, nmax + 1):
        L = 2 * n
        states = core.catalan_states(L)
        S0 = [x for x in states if core.in_S0(x)]
        S1 = [x for x in states if core.in_S1(x)]
        fuss = comb(3 * n, n) // (2 * n + 1)

        enc = {}
        for x in S0:
            b = enc_aligned(x)
            assert is_2dyck(b, n), (x, b)
            assert b not in enc
            enc[b] = x
        bij = len(enc) == fuss == len(S0)

        # A alternates the two sectors; K = A^2 permutes S0
        alt = all(core.in_S1(core.A(x)) for x in S0) and \
              all(core.in_S0(core.A(y)) for y in S1)
        K = {x: core.A(core.A(x)) for x in S0}
        perm = len(set(K.values())) == len(S0) and set(K.values()) == set(S0)

        # no-flat theorem
        flat = sum(1 for x in S0 + S1 if hgt(core.A(x)) == hgt(x))

        # one-bit-carry kernel for the aligned -> staggered half-sweep
        bad = 0
        for x in S0:
            b = enc_aligned(x)
            e, rec = 0, []
            for j in range(n):
                (u, v), e = kappa(e, b[2 * j], b[2 * j + 1])
                rec.append((u, v))
            want = enc_staggered(core.A(x))
            got = [rec[0][0]]
            for j in range(n - 1):
                got.extend([rec[j][1], rec[j + 1][0]])
            got.append(0)
            if tuple(got) != want:
                bad += 1
        print(f"n={n} L={L} |S0|={len(S0):6d} Fuss={fuss:6d} bijection={bij} "
              f"alternates={alt} K_permutation={perm} flat_steps={flat} "
              f"carrier_failures={bad}")


if __name__ == "__main__":
    main()
