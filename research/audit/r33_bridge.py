"""Round-33: the checkpoint reversor Phi = A^(L-2+z) o E, and the reduction.

The round-14 time identity t_k + t_{N-k} = p_L - L + z(x_k) was proved by
telescoping, conditional on the (unproved) checkpoint-gap palindrome.  The
round-26 reversor suggests an unconditional route.  Define, on checkpoints
(states with x_0 = 1, plus E_L with the z(E_L) = L convention),

    Phi(x) = A^(L-2+z(x)) E(x).

Then Phi(E_L) = E_L, and using E A^w = A^-w E plus the ledger fact
z(H(b)) = b_0 (full-state form: z(x_{k+1}) = tau_k, the return time), the
conjugation Phi H = H^-1 Phi algebraically reduces to two pointwise lemmas:

    (A)  Phi maps checkpoints to checkpoints;
    (B)  z(Phi(x)) = z(x).

(A) + (B) + anchor then force Phi(x_k) = x_{N-k}, i.e. the time identity,
with no gap palindrome needed.  This verifies every link at walkable L:

    1. the time identity itself (previously verified only at L = 69);
    2. z(x_{k+1}) = tau_k including the b_0 = 0 edge case;
    3. (A) and (B);
    4. the conjugation Phi H = H^-1 Phi directly.

usage: python3 r33_bridge.py [L ...]
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


def z(x):
    for i, v in enumerate(x):
        if v == 0:
            return i
    return len(x)


def main():
    Ls = [int(a) for a in sys.argv[1:]] or [8, 10, 12, 13, 14, 16, 18, 20, 21, 22]
    print(f"{'L':>3} {'p':>6} {'N':>5}  identity  z=tau   (A)   (B)  conj  anchor")
    for L in Ls:
        st = cycle(L)
        p = len(st)
        # E_L-rooted: st[0] = R shifted... build E_L-rooted cycle
        EL = (1,) * L
        r = st.index(EL)
        st = st[r:] + st[:r]
        idx = {x: t for t, x in enumerate(st)}

        cps = [t for t, x in enumerate(st) if x[0] == 1]
        N = len(cps)
        t_of = cps + [p]           # t_N = p (cyclic return to E_L)

        zs = [z(st[t]) if t > 0 or st[t] != EL else L for t in cps]
        zs[0] = L                  # z(E_L) = L convention

        id_bad = sum(1 for k in range(N)
                     if (t_of[k] + t_of[N - k]) != p - L + zs[k])

        # z(x_{k+1}) = tau_k, full-state form of z(H(b)) = b_0
        ztau_bad = sum(1 for k in range(N - 1)
                       if z(st[t_of[k + 1]]) != t_of[k + 1] - t_of[k])

        # Phi via the cycle: Phi(x_{t}) = x_{sigma(t) + L - 2 + z}
        def phi_time(k):
            return ((p - 2 * L + 2) - t_of[k] + L - 2 + zs[k]) % p

        a_bad = sum(1 for k in range(N) if st[phi_time(k)][0] != 1 and phi_time(k) != 0)
        b_bad = 0
        for k in range(N):
            u = phi_time(k)
            zu = L if u == 0 else z(st[u])
            if zu != zs[k]:
                b_bad += 1

        # conjugation Phi H = H^-1 Phi on checkpoint indices
        pos = {t: j for j, t in enumerate(cps)}
        conj_bad = 0
        for k in range(N):
            u = phi_time(k)
            j = pos.get(u if u != 0 else 0)
            if j is None:
                conj_bad += 1
                continue
            # H advances k -> k+1; want Phi(x_{k+1}) = H^-1 Phi(x_k) = x_{j-1}
            u2 = phi_time((k + 1) % N)
            j2 = pos.get(u2 if u2 != 0 else 0)
            if j2 is None or j2 != (j - 1) % N:
                conj_bad += 1

        anchor = phi_time(0) == 0
        print(f"{L:>3} {p:>6} {N:>5}  {id_bad:>8} {ztau_bad:>6} {a_bad:>5} {b_bad:>5} "
              f"{conj_bad:>5}  {anchor}")


if __name__ == "__main__":
    main()
