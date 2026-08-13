"""Round-26 audit: does the claimed pointed reversor E_n actually exist?

The write-up asserts an explicit recursively computable involution
E_n : C_n -> C_n with

    E_n^2 = id,     E_n A_n = A_n^-1 E_n,     h(E_n x) = h(x),

built from the round-18 ECO lift by

    E_{n+1}(t, r) = ( E_n(A_n t),      r      )   if r <= h(A_n t),
                    ( E_n(Phat_n(t)),  h(t)   )   if r = h(t), h(A_n t) = h(t)-1.

Everything here is rebuilt from `core.A` and the already-verified round-18
lam/unlam/P, so this is an independent check of the construction rather than a
transcription of it.

Note: the mere existence of an involution conjugating A to A^-1 is vacuous
(every finite permutation is a product of two involutions).  The content is
height preservation, which is what forces E_n(R_n) = R_n.

usage: python3 r26_reversor.py [nmax]
"""
import sys

import core
from r18_eco import h, lam, unlam, P


def build(nmax):
    """E_n as an explicit dict, level by level.  Returns {n: {state: state}}."""
    Es = {1: {(1,): (1,)}}
    for n in range(1, nmax):
        prev, cur = Es[n], {}
        for x in core.catalan_states(n + 1):
            t, r = unlam(x)
            At = core.A(t)
            if r <= h(At):
                cur[x] = lam(prev[At], r)
            else:
                # exceptional top fibre: r = h(t), h(At) = h(t) - 1
                assert r == h(t) and h(At) == h(t) - 1, (x, t, r)
                cur[x] = lam(prev[P(t)], h(t))
        Es[n + 1] = cur
    return Es


def audit(n, En):
    states = core.catalan_states(n)
    S = set(states)
    Ainv = {core.A(x): x for x in states}

    closed = sum(1 for x in states if En[x] not in S)
    invol = sum(1 for x in states if En[En[x]] != x)
    conj = sum(1 for x in states if En[core.A(x)] != Ainv[En[x]])
    ht = sum(1 for x in states if h(En[x]) != h(x))

    R = (n,) + (0,) * (n - 1)
    fixesR = En.get(R) == R

    # canonical consequence: E(A^s R) = A^-s R
    walk, bad_walk = R, 0
    back = R
    for _ in range(len(states)):
        walk = core.A(walk)
        back = Ainv[back]
        if En[walk] != back:
            bad_walk += 1
        if walk == R:
            break
    p = 0
    y = R
    while True:
        y = core.A(y)
        p += 1
        if y == R:
            break

    print(f"n={n:2d} |C_n|={len(states):7d} p_n={p:8d} "
          f"not_closed={closed} invol_fail={invol} conj_fail={conj} "
          f"height_fail={ht} E(R)=R:{fixesR} canon_reflect_fail={bad_walk}")
    return closed == invol == conj == ht == bad_walk == 0 and fixesR


def cycle_pairing(n, En):
    """Every A_n-cycle is self-reversing or paired with a reversed partner."""
    states = core.catalan_states(n)
    seen, cycles = set(), []
    for x in states:
        if x in seen:
            continue
        c, y = [], x
        while y not in seen:
            seen.add(y)
            c.append(y)
            y = core.A(y)
        cycles.append(c)
    idx = {}
    for i, c in enumerate(cycles):
        for y in c:
            idx[y] = i
    self_rev = sum(1 for i, c in enumerate(cycles) if idx[En[c[0]]] == i)
    return len(cycles), self_rev, len(cycles) - self_rev


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    Es = build(nmax)
    ok = True
    for n in range(2, nmax + 1):
        ok &= audit(n, Es[n])
    print("\nALL CHECKS PASS" if ok else "\nFAILURES PRESENT")
    for n in (10, 11, 12):
        if n in Es:
            tot, sr, pr = cycle_pairing(n, Es[n])
            print(f"n={n}: {tot} cycles, {sr} self-reversing, {pr} in {pr // 2} reversed pairs")


if __name__ == "__main__":
    main()
