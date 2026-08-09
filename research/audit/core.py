"""Reconstruction of the system from the write-up's own definitions.

State space C_L : x in Z_{>=0}^L, sum(x) = L, and prefix sums  sum_{i<=k} x_i >= k+1.
(This is forced by: E_L = (1,...,1) is a state, sum is L not L-1, and the
 "prefix inequality forces the final coordinate to be at most 1".)

A_L : left-to-right in-place sweep, at each i in 0..L-2 apply the local gate
        x_i = 0        -> nothing
        (1, b)         -> (1+b, 0)
        (a, b), a >= 2 -> (a-1, b+1)
"""
from math import comb


def A(x):
    x = list(x)
    L = len(x)
    for i in range(L - 1):
        a = x[i]
        if a == 0:
            continue
        if a == 1:
            x[i] = 1 + x[i + 1]
            x[i + 1] = 0
        else:
            x[i] = a - 1
            x[i + 1] += 1
    return tuple(x)


def catalan_states(L):
    res = []
    cur = [0] * L

    def rec(pos, s):
        if pos == L:
            if s == L:
                res.append(tuple(cur))
            return
        lo = max(0, pos + 1 - s)
        for v in range(lo, L - s + 1):
            cur[pos] = v
            rec(pos + 1, s + v)
        cur[pos] = 0

    rec(0, 0)
    return res


def E(L):
    return tuple([1] * L)


def is_state(x):
    L = len(x)
    if sum(x) != L:
        return False
    s = 0
    for k, v in enumerate(x):
        s += v
        if s < k + 1:
            return False
    return True


def good(u, v):
    if u == 0:
        return v % 2 == 0
    return u % 2 == 1 and v % 2 == 1


def in_S0(x):
    n = len(x) // 2
    return len(x) % 2 == 0 and all(good(x[2 * j], x[2 * j + 1]) for j in range(n))


def in_S1(x):
    n = len(x) // 2
    if len(x) % 2 or n == 0:
        return False
    if not (x[0] > 0 and x[0] % 2 == 0):
        return False
    if x[2 * n - 1] != 0:
        return False
    return all(good(x[2 * j + 1], x[2 * j + 2]) for j in range(n - 1))


def z(b):
    """leftmost-leaf depth"""
    for k, v in enumerate(b):
        if v == 0:
            return k
    return len(b)


def orbit_period(x0):
    x = x0
    n = 0
    while True:
        x = A(x)
        n += 1
        if x == x0:
            return n


def checkpoints(L):
    """Full canonical A-orbit of E_L, returning (p_L, list of (state, gap))."""
    x0 = E(L)
    x = x0
    p = 0
    cps = []
    last = 0
    while True:
        x = A(x)
        p += 1
        if x[0] == 1:
            cps.append((x, p - last))
            last = p
        if x == x0:
            return p, cps
