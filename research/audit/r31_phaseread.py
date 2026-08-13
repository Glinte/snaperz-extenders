"""Round-31: how many bits of canonical phase are affinely readable from a state?

A cheap phase-reader would be an algorithm: p = t(A^-1 R) + 1.  One bit is
known for even L (the aligned/staggered sector is t mod 2).  This measures how
deep readability goes, using F_q-affine functionals over a cheap feature map as
the proxy class:

    features(x) = [x_i mod q] ++ [x_i == 0] ++ [x_i == 1] ++ [pair parity]
                  ++ [good-pair indicators] ++ [aligned-sector bit] ++ [1]

(all O(L)-computable).  For each prime power q^k | p we ask: is the k-th digit
of t base q an affine functional of the features, *conditionally* on knowing
the lower digits (i.e. restricted to the states where those digits vanish)?
Conditional is the right notion: a recursive reader peels one digit at a time.

The decisive cases are L = 10 (p = 32 = 2^5) and L = 11 (p = 64 = 2^6), where
the phase group is entirely 2-adic: if affine readability stops shallow even
there, phase reading is dead as an attack shape; if it goes to the bottom,
it is the first positive algorithmic signal in the project.

usage: python3 r31_phaseread.py [L ...]
"""
import sys

import core


def cycle(L):
    R = (L,) + (0,) * (L - 1)
    states, x = [], R
    while True:
        states.append(x)
        x = core.A(x)
        if x == R:
            break
    return states


def good_pair(a, b):
    return (a == 0 and b % 2 == 0) or (a % 2 == 1 and b % 2 == 1)


def features(x, q):
    L = len(x)
    f = [v % q for v in x]
    f += [1 if v == 0 else 0 for v in x]
    f += [1 if v == 1 else 0 for v in x]
    if L % 2 == 0:
        pairs = [(x[2 * i], x[2 * i + 1]) for i in range(L // 2)]
        f += [(a + b) % 2 for a, b in pairs]
        f += [1 if good_pair(a, b) else 0 for a, b in pairs]
        f += [1 if all(good_pair(a, b) for a, b in pairs) else 0]
    f += [1]
    return [v % q for v in f]


def solvable(rows, target, q):
    """Is target in the F_q column span of rows?  Gaussian elimination."""
    if not rows:
        return True
    n = len(rows[0])
    aug = [list(r) + [t % q] for r, t in zip(rows, target)]
    inv = [pow(a, q - 2, q) if a else 0 for a in range(q)]
    piv_rows = []
    piv_cols = []
    r = 0
    for c in range(n + 1):
        sel = None
        for i in range(r, len(aug)):
            if aug[i][c]:
                sel = i
                break
        if sel is None:
            continue
        aug[r], aug[sel] = aug[sel], aug[r]
        s = inv[aug[r][c]]
        aug[r] = [(v * s) % q for v in aug[r]]
        for i in range(len(aug)):
            if i != r and aug[i][c]:
                m = aug[i][c]
                aug[i] = [(a - m * b) % q for a, b in zip(aug[i], aug[r])]
        piv_cols.append(c)
        r += 1
        if r == len(aug):
            break
    # solvable iff no pivot in the last (target) column
    return (n not in piv_cols)


def factor(m):
    f, d = {}, 2
    while d * d <= m:
        while m % d == 0:
            f[d] = f.get(d, 0) + 1
            m //= d
        d += 1
    if m > 1:
        f[m] = f.get(m, 0) + 1
    return f


def main():
    Ls = [int(a) for a in sys.argv[1:]] or [7, 10, 11, 12, 13, 14, 16, 18, 20, 21, 22, 24]
    for L in Ls:
        st = cycle(L)
        p = len(st)
        fac = factor(p)
        out = []
        for q, e in sorted(fac.items()):
            if q > 13:
                out.append(f"{q}: (prime too large, skipped)")
                continue
            digits = []
            for k in range(e):
                qk = q ** k
                sub = [t for t in range(p) if t % qk == 0]
                rows = [features(st[t], q) for t in sub]
                tgt = [(t // qk) % q for t in sub]
                ok = solvable(rows, tgt, q)
                digits.append("Y" if ok else "n")
                if not ok:
                    break
            out.append(f"{q}^{e}: digits {' '.join(digits)}")
        print(f"L={L:3d} p={p:6d} = {'*'.join(f'{q}^{e}' if e>1 else str(q) for q,e in sorted(fac.items())):14s} {'  |  '.join(out)}")


if __name__ == "__main__":
    main()
