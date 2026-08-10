"""Round-15 audit, part E: is there a second 16-cell peel?

The R_r renormalization drops the first 16 cells for one macro, worth L/(L-16).
The only way that becomes algorithmically interesting is if the renormalized tail
dynamics re-enters the same shielded form, so that 16, then 32, then 48 cells
could be peeled at successive levels.  That needs the R_r orbit to meet the set
{E_16 || w}.  It does not.
"""
from core import A, E, z
from r15a import R

print("--- does the R_r orbit of E_r stay in the shielded form E_16 || w? ---")
print("   r   first j > 0 with y_j[:16] != E_16   y_1[:16]")
for r in range(17, 31):
    y = E(r)
    bad = None
    for j in range(1, 121):
        y = R(y)
        if y[:16] != (1,) * 16:
            bad = j
            break
    print(f"  {r:2d}   {str(bad):>6}   {R(E(r))[:16]}")

print()
print("--- does it ever come back? ---")
print("   r      M_r   j>0 with E_16 prefix   E_8   E_4   longest run of leading 1s")


def leading_ones(v):
    n = 0
    while n < len(v) and v[n] == 1:
        n += 1
    return n


for r in range(17, 27):
    y0 = E(r)
    y = y0
    seq = []
    while True:
        seq.append(y)
        y = R(y)
        if y == y0:
            break
    tail = seq[1:]
    n16 = sum(1 for v in tail if v[:16] == (1,) * 16)
    n8 = sum(1 for v in tail if v[:8] == (1,) * 8)
    n4 = sum(1 for v in tail if v[:4] == (1,) * 4)
    print(f"  {r:2d}  {len(seq):7d}   {n16:18d}  {n8:4d}  {n4:4d}   "
          f"{max(leading_ones(v) for v in tail):5d}")
