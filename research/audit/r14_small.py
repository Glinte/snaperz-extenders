"""Round-14 audit: the small combinatorial claims.

Uses only core.py's reconstruction of A_L from the gate rule.
"""
import sys
from core import A, E, catalan_states, z, is_state

sys.setrecursionlimit(100000)


# ---------- state <-> plane tree ----------
# x in C_L  <->  plane tree on L+1 vertices, x = preorder degree sequence
# truncated to its first L entries (the last preorder vertex is always a leaf).

def to_tree(x):
    deg = list(x) + [0]
    pos = 0

    def build():
        nonlocal pos
        d = deg[pos]
        pos += 1
        return tuple(build() for _ in range(d))

    t = build()
    assert pos == len(deg)
    return t


def to_state(t):
    out = []

    def walk(u):
        out.append(len(u))
        for c in u:
            walk(c)

    walk(t)
    assert out[-1] == 0
    return tuple(out[:-1])


def theta(t):
    """Dershowitz-Zaks/Deutsch involution: Theta(T1 |x T2) = Theta(T2) |x Theta(T1),
    where T1 |x T2 attaches T2 as the new leftmost child of T1's root."""
    if not t:
        return t
    t2 = t[0]              # first subtree of the root
    t1 = t[1:]             # T1 = t with that subtree removed
    a, b = theta(t2), theta(t1)
    return (b,) + a        # Theta(T2) |x Theta(T1)


def Theta(x):
    return to_state(theta(to_tree(x)))


def H(b):
    """checkpoint map on C_m, via the full state (1,b) in C_{m+1}."""
    x = (1,) + tuple(b)
    for _ in range(b[0] + 1):
        x = A(x)
    assert x[0] == 1, x
    return x[1:]


def cycle(x0):
    orb = [x0]
    x = A(x0)
    while x != x0:
        orb.append(x)
        x = A(x)
    return orb


ok = lambda c: "OK " if c else "FAIL"

# ---------- 1. Theta is an involution and swaps root degree / leftmost depth ----------
for m in (6, 8, 9):
    st = catalan_states(m)
    inv = all(Theta(Theta(x)) == x for x in st)
    swap = all(Theta(x)[0] == z(x) and z(Theta(x)) == x[0] for x in st)
    print(f"[1] m={m:2d} |C_m|={len(st):6d}  Theta involution {ok(inv)}  swaps (b0,z) {ok(swap)}")

# ---------- 2. H = Theta . G with G root-degree preserving; subtree-size multiset ----------
for m in (8, 9, 10):
    st = catalan_states(m)
    zH = all(z(H(b)) == b[0] for b in st)
    # G = Theta . H
    keep = 0
    for b in st:
        g = Theta(H(b))
        assert g[0] == b[0]
        # multiset of root-child subtree sizes (vertex counts)
        a = sorted(len(to_state(c)) + 1 for c in to_tree(b))
        c = sorted(len(to_state(d)) + 1 for d in to_tree(g))
        if a == c:
            keep += 1
    print(f"[2] m={m:2d} |C_m|={len(st):6d}  z(H(b))=b0 {ok(zH)}  "
          f"G preserves root-degree OK  subtree-size multiset preserved: {keep}/{len(st)}")

# ---------- 3. the A_13 canonical cycle and its Theta-compatible phases ----------
E13 = E(13)
orb = cycle(E13)
p13 = len(orb)
idx = {s: i for i, s in enumerate(orb)}
print(f"[3] p_13 = {p13}   (claim 129)  {ok(p13 == 129)}")

compat = [(phi, idx[Theta(y)]) for phi, y in enumerate(orb) if Theta(y) in idx]
law = all(t == (117 - phi) % p13 for phi, t in compat)
print(f"[3] Theta-compatible phases: {len(compat)} of {p13}  (claim 15)  {ok(len(compat) == 15)}")
print(f"    set = {sorted(p for p, _ in compat)}")
print(f"    Theta(y_phi) = y_(117-phi) on all of them {ok(law)}")
survivors = {19, 21, 23, 25}
inter = survivors & {p for p, _ in compat}
print(f"[3] {{19,21,23,25}} cap Theta-compatible = {inter}  (claim {{23}})  {ok(inter == {23})}")

# ---------- 4. y_46 = A^-435 Theta(y_23) ----------
t23 = idx[Theta(orb[23])]
shift = 435 % p13
print(f"[4] Theta(y_23) = y_{t23} (claim 94) {ok(t23 == 94)};  435 mod {p13} = {shift} (claim 48) {ok(shift == 48)}")
print(f"    A^-435 Theta(y_23) = y_{(t23 - shift) % p13} (claim 46) {ok((t23 - shift) % p13 == 46)}")

# ---------- 5. the L=18 checkpoint gap word ----------
x0 = E(18)
x, p, last, gaps = x0, 0, 0, []
while True:
    x = A(x)
    p += 1
    if x[0] == 1:
        gaps.append(p - last)
        last = p
    if x == x0:
        break
print(f"[5] p_18 = {p} (claim 468) {ok(p == 468)};  #gaps = {len(gaps)} (claim 112) {ok(len(gaps) == 112)};  sum = {sum(gaps)}")
print(f"    first 16 gaps = {gaps[:16]}   (claim starts 2,4,2,6,2,4,2,10 repeating)")
import json
import os
json.dump(gaps, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "g18.json"), "w"))
