"""Round-15 audit, part F: the cycle-quotient audit for Gamma_s.

The write-up names this the decisive test for a phase-coordinate algorithm.
Gamma_s = Theta_s . S_{1,s} . A_s is the root-one section of G_{s+1}.  If
Gamma_s carried the canonical A_s-cycle into ONE A_s-cycle with an affine action
on phase, states could be represented as (cycle id, phase mod its period) and the
whole hierarchy would collapse to modular arithmetic.  If instead it scatters one
source cycle across many target cycles, phase compression dies.
"""
from core import A, E, z
from r15a import Gamma


def canonical_cycle(s):
    y0 = E(s)
    orb = [y0]
    y = A(y0)
    while y != y0:
        orb.append(y)
        y = A(y)
    return orb


def cycle_of(x, cache):
    """Identify the A-cycle containing x by its lexicographic minimum."""
    if x in cache:
        return cache[x]
    orb = [x]
    y = A(x)
    while y != x:
        orb.append(y)
        y = A(y)
    key = min(orb)
    for v in orb:
        cache[v] = (key, len(orb))
    return cache[x]


print("--- does Gamma_s map the canonical A_s-cycle into one A_s-cycle? ---")
print("   s    p_s   targets  per point   biggest   shortest  longest   phase map affine?")
for s in range(2, 18):
    orb = canonical_cycle(s)
    ps = len(orb)
    idx = {y: i for i, y in enumerate(orb)}
    cache = {}
    targets = {}
    onto_self = 0
    for y in orb:
        g = Gamma(y)
        key, ln = cycle_of(g, cache)
        targets.setdefault((key, ln), 0)
        targets[(key, ln)] += 1
        if g in idx:
            onto_self += 1

    # if it lands entirely on the canonical cycle, is phi -> a*phi + b affine?
    affine = "n/a"
    if onto_self == ps:
        img = [idx[Gamma(y)] for y in orb]
        found = None
        for a in range(ps):
            b = (img[0] - a * 0) % ps
            if all(img[t] == (a * t + b) % ps for t in range(ps)):
                found = (a, b)
                break
        affine = f"yes {found}" if found else "NO"

    lens = [ln for (_, ln) in targets]
    print(f"  {s:2d}  {ps:5d}   {len(targets):6d}   {len(targets) / ps:8.2f}   "
          f"{max(targets.values()):6d}   {min(lens):6d} {max(lens):8d}   {affine}")
