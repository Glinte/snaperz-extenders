"""Round-18/19 audit: the rightmost-peak (ECO) lift, the fibre transition
theorem, the return permutation Xi, and the two-slot monodromy.

usage: python3 r18_eco.py [nmax] [Lmax]
"""
import sys

import core


def last_pos(t):
    m = len(t) - 1
    while m >= 0 and t[m] == 0:
        m -= 1
    return m


def h(t):
    return len(t) - last_pos(t)


def lam(t, r):
    """Insert one peak into the final descent, leaving r old down-steps after."""
    n = len(t)
    m = last_pos(t)
    hh = n - m
    assert 0 <= r <= hh
    if r == hh:
        return tuple(t[:m]) + (t[m] + 1,) + (0,) * hh
    return tuple(t[:m + 1]) + (0,) * (hh - r - 1) + (1,) + (0,) * r


def unlam(s):
    """Inverse of lam: delete the rightmost peak.  Returns (t, r)."""
    n = len(s)
    m = last_pos(s)
    if s[m] == 1:
        # the rightmost peak is its own coordinate
        r = n - 1 - m
        t = tuple(s[:m]) + (0,) * (n - 1 - m)
        t = t[:n - 1]
        # strip: coordinate m disappears
        t = tuple(s[:m]) + (0,) * (n - 1 - m)
        return t, r
    t = tuple(s[:m]) + (s[m] - 1,) + (0,) * (n - 1 - m - 1)
    return t, h(t)


def P(t):
    """Gates strictly before the last positive coordinate."""
    x = list(t)
    m = last_pos(t)
    for i in range(m):
        a = x[i]
        if a == 1:
            x[i] = 1 + x[i + 1]
            x[i + 1] = 0
        elif a >= 2:
            x[i] = a - 1
            x[i + 1] += 1
    return tuple(x)


def check_lift(n):
    """lam is a bijection, and the fibre transition theorem holds."""
    states = core.catalan_states(n)
    up = set(core.catalan_states(n + 1))
    seen = {}
    for t in states:
        for r in range(h(t) + 1):
            s = lam(t, r)
            assert s in up, (t, r, s)
            assert s not in seen
            seen[s] = (t, r)
            assert unlam(s) == (t, r), (t, r, s, unlam(s))
    assert len(seen) == len(up), (len(seen), len(up))

    census = {-1: 0, 0: 0, 1: 0}
    bad = 0
    exc = 0
    for t in states:
        At = core.A(t)
        d = h(At) - h(t)
        assert d in (-1, 0, 1), (t, d)
        census[d] += 1
        hh = h(t)
        for r in range(hh + 1):
            got = unlam(core.A(lam(t, r)))
            if d == 1:
                want = (At, r) if r < hh else (At, hh + 1)
            elif d == 0:
                want = (At, r) if r < hh - 1 else (At, hh if r == hh - 1 else hh - 1)
            else:
                if r < hh - 2:
                    want = (At, r)
                elif r == hh - 2:
                    want = (At, hh - 1)
                elif r == hh - 1:
                    want = (At, hh - 2)
                else:
                    want = (P(t), hh - 1)
                    exc += 1
            if got != want:
                bad += 1
    cn1 = core.comb(2 * (n - 1), n - 1) // n if n >= 2 else 0
    print(f"n={n:2d} |C_n|={len(states):6d} lift_ok fibre_failures={bad} "
          f"up={census[1]} down={census[-1]} (C_(n-1)={cn1}) flat={census[0]} "
          f"exceptional_edges={exc}")
    return bad == 0 and census[1] == census[-1] == cn1


def eta(v):
    return lam(lam(v, h(v)), h(v))


def xi(m):
    """Return map on the section eta(C_m) inside C_{m+2}."""
    sec = {}
    for v in core.catalan_states(m):
        sec[eta(v)] = v
    succ, wt = {}, {}
    for s, v in sec.items():
        x, w = core.A(s), 1
        while x not in sec:
            x = core.A(x)
            w += 1
        succ[v] = sec[x]
        wt[v] = w
    assert len(set(succ.values())) == len(succ), "Xi is not a permutation"
    return succ, wt


def canonical(m):
    succ, wt = xi(m)
    R = (m,) + (0,) * (m - 1)
    assert eta(R) == core.A((m + 2,) + (0,) * (m + 1)), "anchor mismatch"
    word, v = [], R
    while True:
        word.append(wt[v])
        v = succ[v]
        if v == R:
            break
    return word


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    Lmax = int(sys.argv[2]) if len(sys.argv) > 2 else 26
    for n in range(2, nmax + 1):
        assert check_lift(n)

    print()
    for m in range(1, min(nmax, 8) + 1):
        word = canonical(m)
        L = m + 2
        p = sum(word)
        pal = all(word[i] == word[len(word) - 1 - i] for i in range(len(word)))
        print(f"L={L:2d} Xi-period={len(word):6d} sum(w)={p:8d} p_L={core.orbit_period(core.E(L)):8d} "
              f"palindrome={pal}")

    print("\ncanonical Xi weight words:")
    for m in (3, 5, 11):
        if m + 2 <= Lmax:
            print(f"  L={m+2}: {canonical(m)}")


if __name__ == "__main__":
    main()
