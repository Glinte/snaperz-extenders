from core import *
import sys

def H(x):
    y = A(x); t = 1
    while y[0] != 1:
        y = A(y); t += 1
    return y, t

def Hb(b):
    """H on the reduced state b (checkpoint x = (1,b))."""
    y, t = H((1,) + tuple(b))
    return y[1:], t

# ---------- trees ----------
def to_tree(b):
    deg = list(b) + [0]
    it = iter(range(len(deg)))
    pos = [0]
    def build():
        i = pos[0]; pos[0] += 1
        return tuple(build() for _ in range(deg[i]))
    t = build()
    assert pos[0] == len(deg)
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

def Theta(t):
    if not t:
        return ()
    c1, rest = t[0], t[1:]
    a = Theta(c1)          # Theta(T2)
    b = Theta(rest)        # Theta(T1)
    return (b,) + a        # Theta(T2) join Theta(T1): b becomes leftmost child of a's root

def ThetaS(b):
    return to_state(Theta(to_tree(b)))

print('--- SS5: Theta involution and root/depth exchange ---')
for m in range(1, 10):
    S = catalan_states(m)
    inv = all(ThetaS(ThetaS(b)) == b for b in S)
    exch = all(ThetaS(b)[0] == z(b) and z(ThetaS(b)) == b[0] for b in S)
    print(f'm={m:>2} |C_m|={len(S):>5} involution:{inv} (Theta b)_0=z(b) & z(Theta b)=b_0: {exch}')

print()
print('--- SS5: z(H_m(b)) = b_0  over the ENTIRE state space ---')
for m in range(1, 10):
    S = [b for b in catalan_states(m)]
    ok = all(z(Hb(b)[0]) == b[0] for b in S)
    # is H a permutation of C_m?  and does Theta H Theta = H^{-1}?
    img = {Hb(b)[0] for b in S}
    perm = len(img) == len(S)
    dual = all(ThetaS(Hb(ThetaS(b))[0]) == b for b in S) if perm else None
    # G = Theta o H preserves root degree
    grd = all(ThetaS(Hb(b)[0])[0] == b[0] for b in S)
    print(f'm={m:>2} z(H b)=b_0:{ok}  H perm:{perm}  (Theta H Theta)=H^-1:{dual}  G preserves root deg:{grd}')

print()
print('--- SS6: the 16-clock  H^112 (E_15, y) = (E_15, R(y)) ---')
for r in range(1, 7):
    Cr = catalan_states(r)
    gapsets = None
    ok = True
    info = []
    for y in Cr:
        b = tuple([1]*15) + tuple(y)
        gaps = []
        cur = b
        for _ in range(112):
            cur, g = Hb(cur)
            gaps.append(g)
        if cur[:15] != tuple([1]*15):
            ok = False
        yp = cur[15:]
        # A_r^435 y
        w = y
        for _ in range(435 % (r and 1 or 1) * 0 + 435):
            w = A(w)
        c = w[0]
        info.append((sum(gaps[:111]), gaps[111], c, z(yp), sum(gaps)))
        if gapsets is None:
            gapsets = gaps[:111]
        elif gapsets != gaps[:111]:
            ok = False
    s450 = {i[0] for i in info}
    last_ok = all(i[1] == 16 + i[2] for i in info)
    tot_ok = all(i[4] == 466 + i[2] for i in info)
    zid = all(i[3] == i[2] for i in info)
    print(f'r={r:>2} |C_r|={len(Cr):>4} prefix restored:{ok}  sum(first 111)={s450}  '
          f'last gap=16+c:{last_ok}  total=466+c:{tot_ok}  z(R(y))=(A^435 y)_0:{zid}')
    sys.stdout.flush()
