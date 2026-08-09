from core import *
from math import comb
import sys

def gen_S0(n):
    """Directly enumerate good-pair states of size 2n (sum 2n, Catalan prefix)."""
    L = 2 * n
    out = []
    cur = [0] * L

    def rec(j, s):
        if j == n:
            if s == L:
                out.append(tuple(cur))
            return
        rem = L - s
        # pair (u,v) good, prefix conditions at 2j and 2j+1
        for u in range(0, rem + 1):
            if s + u < 2 * j + 1:
                continue
            if u != 0 and u % 2 == 0:
                continue
            vs = range(0, rem - u + 1, 2) if u == 0 else range(1, rem - u + 1, 2)
            for v in vs:
                if s + u + v < 2 * j + 2:
                    continue
                cur[2 * j] = u
                cur[2 * j + 1] = v
                rec(j + 1, s + u + v)
        cur[2 * j] = cur[2 * j + 1] = 0

    rec(0, 0)
    return out


def H(x):
    """return map to the next checkpoint (x_0 == 1); returns (image, gap)"""
    y = A(x)
    t = 1
    while y[0] != 1:
        y = A(y)
        t += 1
    return y, t


def cycle_structure(states):
    idx = {s: i for i, s in enumerate(states)}
    nxt = [idx[H(s)[0]] for s in states]
    seen = [False] * len(states)
    cycles = []
    for i in range(len(states)):
        if seen[i]:
            continue
        path = []
        j = i
        while not seen[j]:
            seen[j] = True
            path.append(j)
            j = nxt[j]
        if j in path:
            cycles.append(len(path) - path.index(j))
    return cycles


print(f"{'n':>2} {'#cps':>7} {'N_n':>7} {'canon':>6} {'maxcyc':>7}  tau=x1+1  perm")
for n in range(1, 10):
    S0 = gen_S0(n)
    assert len(S0) == comb(3 * n, n) // (2 * n + 1)
    cps = [x for x in S0 if x[0] == 1]
    Nn = 2 * comb(3 * n - 1, n - 1) // (3 * n - 1)
    imgs = [H(x) for x in cps]
    perm = len({i[0] for i in imgs}) == len(cps) and all(i[0][0] == 1 for i in imgs)
    tau_ok = all(g == x[1] + 1 for x, (_, g) in zip(cps, imgs))
    # canonical H-period from E_{2n}
    x = E(2 * n)
    per = 0
    while True:
        x = H(x)[0]
        per += 1
        if x == E(2 * n):
            break
    cyc = cycle_structure(cps) if len(cps) < 300000 else []
    print(f"{n:>2} {len(cps):>7} {Nn:>7} {per:>6} {max(cyc):>7}  {tau_ok}      {perm}")
    sys.stdout.flush()
