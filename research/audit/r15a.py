"""Round-15 audit: the state-level renormalisation R_r = H . S_30 . A^435.

Everything is rebuilt from core.A (the gate rule) and the tree involution; no
claim from the write-up is used as an input.  Run from inside audit/.
"""
import sys
from core import A, E, catalan_states, z, is_state
from r14_small import Theta, to_tree, to_state

sys.setrecursionlimit(100000)


def ok(b):
    return "OK" if b else "FAIL"


# ---------- the maps ----------

def Araw(x, n=1):
    """The sweep as a raw gate pass (the write-up applies it to non-states)."""
    for _ in range(n):
        x = A(x)
    return x


def J(u):
    """J_m(u) = A_m(u + e_0)."""
    u = (u[0] + 1,) + tuple(u[1:])
    return A(u)


def Jp(u, c):
    for _ in range(c):
        u = J(u)
    return u


def H_sim(b):
    """Checkpoint map, by definition: first return of (1,b) to a_0 == 1."""
    x = (1,) + tuple(b)
    x = A(x)
    t = 1
    while x[0] != 1:
        x = A(x)
        t += 1
    return x[1:], t


def H(b):
    """Closed form H_m(b) = J_m^{b_0}(0, A_{m-1}(b_1..))."""
    return Jp((0,) + Araw(tuple(b[1:])), b[0])


def S(w, t):
    """S_{t,r}(c,v) = (c, A_{r-1}^t v): the forest clock, root held fixed."""
    return (w[0],) + Araw(tuple(w[1:]), t)


def R(y):
    """R_r = H_r . S_{30,r} . A_r^435."""
    return H(S(Araw(tuple(y), 435), 30))


def R_naive(y):
    """The tempting simplification H_r . A_r^435 (no forest phase debt)."""
    return H(Araw(tuple(y), 435))


def R_explicit(y):
    """R_r(y) = J_r^c(0, A_{r-1}^31 v)  with  (c,v) = A_r^435(y)."""
    w = Araw(tuple(y), 435)
    c, v = w[0], w[1:]
    return Jp((0,) + Araw(v, 31), c)


# ---------- 0. the closed form for H, and R's two spellings ----------
print("--- 0. H closed form and the two spellings of R ---")
for m in range(1, 10):
    st = catalan_states(m)
    a = all(H(b) == H_sim(b)[0] for b in st)
    b_ = all(H_sim(b)[1] == b[0] + 1 for b in st)
    print(f"  m={m:2d} |C_m|={len(st):5d}  H=J^b0(0,A b1..) {ok(a)}   return time b_0+1 {ok(b_)}")

same = all(R(y) == R_explicit(y) for r in range(1, 9) for y in catalan_states(r))
print(f"  R = H.S_30.A^435  ==  J^c(0,A^31 v) on all C_r, r<=8: {ok(same)}")
print()


# ---------- 1. R against the literal 112 checkpoint returns ----------
print("--- 1. H_{15+r}^112(E_15 || y) = E_15 || R_r(y) ---")
for r in range(1, 8):
    st = catalan_states(r)
    bad = []
    for y in st:
        b = (1,) * 15 + tuple(y)
        assert is_state((1,) + b)
        x = b
        for _ in range(112):
            x = H_sim(x)[0]
        if x != (1,) * 15 + R(y):
            bad.append(y)
    print(f"  r={r} |C_r|={len(st):4d}  mismatches {len(bad)} {ok(not bad)}")
print()


# ---------- 2. the full-sweep statement and the elapsed time 466+c ----------
print("--- 2. A_{16+r}^{466+c}(E_16 || y) = E_16 || R_r(y),  c = (A_r^435 y)_0 ---")
for r in range(1, 7):
    st = catalan_states(r)
    bad = tbad = 0
    for y in st:
        c = Araw(tuple(y), 435)[0]
        x = (1,) * 16 + tuple(y)
        x = Araw(x, 466 + c)
        if x != (1,) * 16 + R(y):
            bad += 1
        # is 466+c really the *first* return to the shielded prefix?
        u = (1,) * 16 + tuple(y)
        for t in range(1, 466 + c):
            u = A(u)
            if u[:16] == (1,) * 16:
                tbad += 1
                break
    print(f"  r={r} |C_r|={len(st):4d}  state mismatches {bad} {ok(not bad)}   "
          f"earlier prefix returns {tbad}")
print()


# ---------- 3. the naive simplification, and where it first fails ----------
print("--- 3. R =? H . A^435  (the omitted forest clock) ---")
for r in range(1, 9):
    st = catalan_states(r)
    bad = [y for y in st if R(y) != R_naive(y)]
    print(f"  r={r} |C_r|={len(st):5d}  disagreements {len(bad):5d} "
          f"{'(agrees)' if not bad else '<- first failure' if r == 6 else ''}")
print()


# ---------- 4. periods M_r and the reduced recurrence for p ----------
print("--- 4. p_{16+r} = 466 M_r + sum_j z(y_j),  y_{j+1} = R_r(y_j), y_0 = E_r ---")


def orbit_M(r):
    y0 = E(r)
    y = y0
    seq = []
    while True:
        seq.append(y)
        y = R(y)
        if y == y0:
            return len(seq), seq


def p_direct(L):
    x0 = E(L)
    x = A(x0)
    p = 1
    while x != x0:
        x = A(x)
        p += 1
    return p


claimed_M = {1: 1, 2: 1, 3: 2, 4: 1, 5: 4, 6: 2, 7: 4, 8: 2, 9: 8, 10: 10,
             11: 12, 12: 16, 13: 24, 14: 10, 16: 14, 18: 19}
print("   r    L      M_r  claim   466M+sum z   p_L direct   B_L = p-L+1")
for r in list(range(1, 19)):
    L = 16 + r
    M, seq = orbit_M(r)
    pred = 466 * M + sum(z(y) for y in seq)
    # cross-check with the raw sum of macro lengths 466+c_j
    pred2 = sum(466 + Araw(y, 435)[0] for y in seq)
    direct = p_direct(L) if L <= 34 else None
    cm = claimed_M.get(r, "-")
    flag = "" if direct is None else ok(direct == pred)
    assert pred == pred2, (r, pred, pred2)
    print(f"  {r:2d}   {L:2d}   {M:6d}  {str(cm):>5}   {pred:10d}   "
          f"{str(direct):>10}   {pred - L + 1:10d}  {flag}")
print()


# ---------- 5. the reflection candidate Psi = H . Theta ----------
print("--- 5. Psi_m = H_m . Theta_m as the checkpoint reflection ---")


def Rstate(L):
    """The retraction state: A_L^{L-1}(R_L) = E_L."""
    # walk the canonical cycle backwards L-1 steps from E_L
    x = E(L)
    p = p_direct(L)
    return Araw(x, (p - (L - 1)) % p)


for m in range(1, 8):
    Psi = lambda x: H(Theta(x))
    fixE = Psi(E(m)) == E(m)
    st = catalan_states(m)
    keepz = all(z(Psi(x)) == z(x) for x in st)
    # does it reflect the canonical H-cycle:  Psi(H^k E) = H^{-k} E ?
    orb = [E(m)]
    while True:
        nxt = H(orb[-1])
        if nxt == E(m):
            break
        orb.append(nxt)
    n = len(orb)
    refl = all(Psi(orb[k]) == orb[(-k) % n] for k in range(n))
    print(f"  m={m} |cycle|={n:3d}  Psi(E)=E {ok(fixE)}  z-preserving {ok(keepz)}  "
          f"reflects cycle {ok(refl)}")

b = H(E(6))
print(f"  m=6 witness:  H_6(E_6) = {b}")
print(f"                H_6 Theta_6 (b) = {H(Theta(b))}")
print(f"                R_6 = {Rstate(6)}   (equal: {H(Theta(b)) == Rstate(6)})")
print()


# ---------- 6. the small-tail phase map Gamma_s ----------
print("--- 6. Gamma_s = Theta_s . S_{1,s} . A_s  as the root-one section of G ---")


def G(x):
    return Theta(H(x))


def Gamma(v):
    """claimed: Gamma_s = Theta_s . S_{1,s} . A_s   (acting on v of length s)"""
    return Theta(S(A(tuple(v)), 1))


for s in range(1, 9):
    # root-one section: G_{s+1}(1,v) = (1, Gamma_s(v)) for v with (1,v) in C_{s+1}
    vs = [x[1:] for x in catalan_states(s + 1) if x[0] == 1]
    sec = all(G((1,) + v)[0] == 1 for v in vs)
    match = all(G((1,) + v)[1:] == Gamma(v) for v in vs)
    print(f"  s={s} |root-one fibre|={len(vs):5d}  G keeps root 1 {ok(sec)}  "
          f"section = Theta.S_1.A {ok(match)}")
