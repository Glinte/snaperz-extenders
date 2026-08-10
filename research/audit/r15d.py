"""Round-15 audit, part D: is the new p-recurrence the logged 112-template?

The ledger already has the 112-block Toeplitz template: for m=15..30 the
checkpoint gaps of the canonical orbit fall into blocks of 112 whose first 111
gaps are universal (independent of L) and whose 112th gap ("the hole") varies.
The template "gives weights, never states".

The round-15 write-up gives  p_{16+r} = 466 M_r + sum_j z(y_j)  with y_{j+1} =
R_r(y_j).  If the two are the same object then the block count must be M_r and
the hole must be a function of y_j.  That is what this checks.
"""
from collections import Counter

from core import A, E, z
from r15a import R, orbit_M


def gaps(L):
    """Checkpoint gaps around the canonical orbit of E_L, starting at E_L."""
    x0 = E(L)
    x = x0
    out = []
    last = 0
    p = 0
    while True:
        x = A(x)
        p += 1
        if x[0] == 1:
            out.append(p - last)
            last = p
        if x == x0:
            return out


print("--- the 112-block template vs the R_r orbit ---")
print("   L   r   #gaps  blocks  M_r   universal 111-prefix?  sum(111)  holes = 16+z(y_j)?")
ref = None
for r in range(1, 19):
    L = 16 + r
    g = gaps(L)
    assert len(g) % 112 == 0, (L, len(g))
    nb = len(g) // 112
    M, seq = orbit_M(r)
    blocks = [g[112 * j:112 * (j + 1)] for j in range(nb)]
    heads = {tuple(b[:111]) for b in blocks}
    uni = len(heads) == 1
    if ref is None:
        ref = next(iter(heads))
    same_as_ref = heads == {ref}
    holes = [b[111] for b in blocks]
    # the write-up's y_j are read at the E_16-prefix returns, i.e. block starts
    pred = [16 + z(y) for y in seq]
    # align: the block sequence must be the z-sequence up to rotation
    ok_hole = nb == M and any(holes == pred[i:] + pred[:i] for i in range(M))
    print(f"  {L:2d}  {r:2d}   {len(g):5d}   {nb:5d}  {M:4d}   "
          f"{'yes' if uni and same_as_ref else 'NO':>19}   {sum(ref):7d}   "
          f"{'yes' if ok_hole else 'NO'}")

print()
print(f"universal 111-gap prefix, distribution: {dict(sorted(Counter(ref).items()))}")
print(f"  count {len(ref)}, sum {sum(ref)}")
