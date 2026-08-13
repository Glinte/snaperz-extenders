"""Round-32: linear complexity of the canonical orbit over odd characteristic.

The ledger's LC = p_L result (F_2, rounds 16-23) kills polynomial-dimensional
linearization over F_2 only.  A d-dimensional linear representation over Q --
the last measurable "cheap composition" shape -- would force every integer
observable of the orbit to satisfy an order-d recurrence over Q, hence over
F_q for all but finitely many primes q.  So Berlekamp-Massey on the integer
sequence (x_t)_0 mod q, for a few odd and large primes, decides it: LC ~ p at
one large prime rules out any small linear representation over characteristic
zero as well.

usage: python3 r32_lcq.py [L ...]
"""
import sys

import core


def cycle_coord0(L):
    R = (L,) + (0,) * (L - 1)
    seq, x = [], R
    while True:
        seq.append(x[0])
        x = core.A(x)
        if x == R:
            break
    return seq


def berlekamp_massey(s, q):
    """Linear complexity of s over F_q."""
    n = len(s)
    c = [1] + [0] * n
    b = [1] + [0] * n
    lc, m, bd = 0, 1, 1
    inv = [0] * q
    for a in range(1, q):
        inv[a] = pow(a, q - 2, q)
    for i in range(n):
        d = s[i] % q
        for j in range(1, lc + 1):
            d = (d + c[j] * s[i - j]) % q
        if d == 0:
            m += 1
        elif 2 * lc <= i:
            t = c[:]
            coef = (d * inv[bd]) % q
            for j in range(0, n - m + 1):
                c[j + m] = (c[j + m] - coef * b[j]) % q
            lc, b, bd, m = i + 1 - lc, t, d, 1
        else:
            coef = (d * inv[bd]) % q
            for j in range(0, n - m + 1):
                c[j + m] = (c[j + m] - coef * b[j]) % q
            m += 1
    return lc


def main():
    Ls = [int(a) for a in sys.argv[1:]] or [13, 14, 16, 18, 20, 21, 22]
    primes = [2, 3, 5, 101, 997]
    print(f"{'L':>3} {'p':>6}  " + "  ".join(f"LC_mod_{q:<3}" for q in primes))
    for L in Ls:
        seq = cycle_coord0(L)
        p = len(seq)
        s2 = seq + seq  # two periods so BM sees the full recurrence
        row = [berlekamp_massey(s2, q) for q in primes]
        print(f"{L:>3} {p:>6}  " + "  ".join(f"{lc:>10}" for lc in row))


if __name__ == "__main__":
    main()
