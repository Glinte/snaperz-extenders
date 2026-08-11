"""Round-22 audit: is the sweep really the claimed plane-tree recursion?

Checks, all against `core.A` rebuilt from the gate rule:
  1. F(T) recursion == A on every Catalan state through n = 12;
  2. the recursive inverse F^-1;
  3. the subtree-endpoint rotor formula for a single gate;
  4. ord(T_i) = lcm(1, ..., n-i);
  5. the append/unary-spine interface law F(U (+) V);
  6. the persistent proper-subtree cache counts D_n.

usage: python3 r22_tree.py [nmax]
"""
import sys
from math import lcm

import core

# ---------------------------------------------------------------- trees


def decode(x):
    """Catalan state -> nested tuple tree (preorder outdegree sequence + 0)."""
    seq = list(x) + [0]
    pos = 0

    def rec():
        nonlocal pos
        d = seq[pos]
        pos += 1
        return tuple(rec() for _ in range(d))

    t = rec()
    assert pos == len(seq)
    return t


def encode(t):
    out = []

    def rec(u):
        out.append(len(u))
        for c in u:
            rec(c)

    rec(t)
    assert out[-1] == 0
    return tuple(out[:-1])


def F(t):
    if not t:
        return ()
    if len(t) == 1:
        raw = ((),) + t[0]
    else:
        raw = (t[0] + (t[1],),) + t[2:]
    return tuple(F(c) for c in raw)


def Finv(y):
    if not y:
        return ()
    z = tuple(Finv(c) for c in y)
    if z[0] == ():
        return (z[1:],)
    return (z[0][:-1],) + (z[0][-1],) + z[1:]


# ---------------------------------------------------------------- gates


def gate(x, i):
    """One local gate at position i (the sweep is gate 0, 1, ..., n-2)."""
    x = list(x)
    a = x[i]
    if a == 1:
        x[i] = 1 + x[i + 1]
        x[i + 1] = 0
    elif a >= 2:
        x[i] = a - 1
        x[i + 1] += 1
    return tuple(x)


def endpoints(x):
    """e_j = last preorder index of the subtree rooted at vertex j."""
    seq = list(x) + [0]
    n = len(seq)
    e = [0] * n
    stack = []
    for j in range(n):
        stack.append([j, seq[j]])
        while stack and stack[-1][1] == 0:
            v = stack.pop()[0]
            e[v] = j
            if stack:
                stack[-1][1] -= 1
    return e


def rotor(e, i):
    """The claimed endpoint update: only e_{i+1} changes."""
    e = list(e)
    if e[i] == i:
        return e
    if e[i + 1] < e[i]:
        e[i + 1] = e[e[i + 1] + 1]
    else:
        e[i + 1] = i + 1
    return e


# ---------------------------------------------------------------- checks


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    total = 0
    for n in range(1, nmax + 1):
        states = core.catalan_states(n)
        total += len(states)
        bad_f = bad_inv = bad_e = 0
        for x in states:
            t = decode(x)
            y = F(t)
            if encode(y) != core.A(x):
                bad_f += 1
            if Finv(y) != t or F(Finv(t)) != t:
                bad_inv += 1
            if n <= 9:
                e = endpoints(x)
                for i in range(n - 1):
                    if rotor(e, i) != endpoints(gate(x, i)):
                        bad_e += 1
        print(f"n={n:2d} states={len(states):7d} sweep_fail={bad_f} "
              f"inverse_fail={bad_inv} rotor_fail={bad_e}")
        assert not (bad_f or bad_inv or bad_e)
    print(f"total Catalan states checked: {total}")

    # gate orders
    print("\nord(T_i) vs lcm(1..n-i):")
    for n in range(2, 9):
        states = core.catalan_states(n)
        idx = {x: k for k, x in enumerate(states)}
        for i in range(n - 1):
            seen = [False] * len(states)
            lens = set()
            for s, x in enumerate(states):
                if seen[s]:
                    continue
                cyc, y = 0, x
                while True:
                    seen[idx[y]] = True
                    y = gate(y, i)
                    cyc += 1
                    if y == x:
                        break
                lens.add(cyc)
            got, want = lcm(*lens) if len(lens) > 1 else lens.pop(), lcm(*range(1, n - i + 1))
            flag = "ok" if got == want else "MISMATCH"
            print(f"  n={n} i={i}: ord={got} lcm(1..{n-i})={want} {flag}")

    # append law
    print("\nappend law F(U (+) V):")
    pairs = bad = 0
    trees = {k: [decode(x) for x in core.catalan_states(k)] for k in range(1, 10)}
    trees[0] = [()]
    for su in range(0, 10):
        for sv in range(1, 10 - su):
            for u in trees[su]:
                for v in trees[sv]:
                    # U (+) V: append V as the last child of the root of U
                    joined = u + (v,)
                    pairs += 1
                    lhs = F(joined)
                    if len(u) >= 2:
                        rhs = F(u) + (F(v),)
                    elif len(u) == 1:
                        rhs = (F(u[0] + (v,)),)
                    else:
                        rhs = ((),) + tuple(F(c) for c in v)
                    if lhs != rhs:
                        bad += 1
    print(f"  pairs={pairs} failures={bad}")

    # persistent proper-subtree cache
    print("\npersistent proper-subtree transforms:")
    for n in (13, 16, 21, 25):
        cache = {}
        x = core.E(n)
        x0, p, new = x, 0, 0

        def Fc(t):
            nonlocal new
            if not t:
                return ()
            if len(t) == 1:
                raw = ((),) + t[0]
            else:
                raw = (t[0] + (t[1],),) + t[2:]
            out = []
            for c in raw:
                if c in cache:
                    out.append(cache[c])
                else:
                    r = Fc(c)
                    cache[c] = r
                    new += 1
                    out.append(r)
            return tuple(out)

        while True:
            t = decode(x)
            y = Fc(t)
            x = encode(y)
            p += 1
            if x == x0:
                break
        print(f"  n={n:2d} p={p:6d} flat_gates={(n-1)*p:9d} new_transforms={new:8d} "
              f"per_sweep={new/p:.3f} gates_per_transform={(n-1)*p/new:.2f}")


if __name__ == "__main__":
    main()
