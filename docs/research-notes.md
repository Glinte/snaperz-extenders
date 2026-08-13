# Exact pulse counts for P = 12 Snaperz extenders

Notes from an attempt to compute `T_L`, the number of pulses a period-12 Snaperz
extender of length `L` takes to complete a cycle, without simulating all of it.
The attempt did not succeed: no subexponential algorithm was found, and the
delivered speedups are constant factors. What follows is the theory that *is*
settled, the data, and — at least as valuable — the routes that were measured
dead, so nobody has to pay for them twice.

Everything here is reproducible from the programs in [`../research`](../research);
that directory's README says which program checks which claim.

## Settled theory

A state is the tuple of air-separated segment lengths `(a_0, ..., a_{n-1})` with
`sum(a) == n`. The reachable states are exactly the **Catalan states**
`a_0 + ... + a_k >= k + 1`, i.e. Dyck paths of semilength `n`.

**Blockless sweep.** Deleting the extender block leaves a sweep `A_L` that runs
left to right over the Catalan states with the local gate

    (0, b)          -> no-op
    (1, b)          -> (1 + b, 0)
    (a, b), a >= 2  -> (a - 1, b + 1)

`A_L` is a permutation. Write `E_L = (1, ..., 1)` for extended and
`R_L = (L, 0, ..., 0)` for retracted; they lie on one `A_L`-cycle of length
`p_L`. Since `A_L^{L-1}(R_L) = E_L`,

    p_L = B_L + (L - 1)

where `B_L` is the first-hitting time `E_L -> R_L`.

**Good-pair / fibre theorem.** Re-inserting the block adds one fibre coordinate
`r`. Tracking the real block together with an auxiliary block, exactly one of
them retracts on the first trip, and

    T_L = (1 + eps_L) * p_L - (L - 1),   so   T_L in {B_L, 2 B_L + L - 1}.

`eps_L = 1` for L = 2, 4, 6, 12, 20, 21, 41 (7 of the first 46). This doubles as
the termination proof. It is paper-level, not formalized.

**Physical -> blockless reduction** is the last-peak lemma
`pi(A_{n+1}(lam(t, r))) = A_n(t)`, verified to n = 6 — but only for `r < h(t)`.
The exceptional set at `r = h(t)` is exactly Catalan(n-1) states, so this must
not be restated as an unconditional `delta . F_L = A_L . delta`.

**Checkpoints** are the states with `a_0 == 1`; `H` is the induced map on the
tail, with return time `b_0 + 1`. Exactly

    B_L = 1 + sum_{j < Q_L} (b^(j)_0 + 1),

with `b^(0) = E_{L-1}` and `Q_L = min{ q : H^q(E_{L-1}) = R_{L-1} }` (checked
L = 3..20). Closed form: `H_m(b) = J_m^{b_0}(0, A_{m-1}(b_1..))` with
`J_m(u) = A_m(u + e_0)`.

**`z(H_m(b)) = b_0`**, where `z(b) = min{ i : b_i = 0 }` — zero counterexamples
over all 23,713 Catalan states to m = 10. Equivalently `(H^{-1}(b))_0 = z(b)`;
in plane-tree terms, leftmost-leaf depth after `H` equals root degree before.
This is the cleanest new structural fact the project produced. It also gives
`H_m = Theta . G_m`, where `Theta` is the plane-tree involution swapping root
degree and leftmost-leaf depth and `G_m` preserves root degree (verified over
all of `C_m`, m <= 9). Note `Theta H Theta = H^{-1}` is **false**.

**Cumulative coordinates and a 3-state carry automaton.** On
`s_i = sum_{j <= i} a_j` each gate writes one cell: `a_i = 0` is a no-op,
`a_i = 1` sets `s_i = s_{i+1}`, `a_i >= 2` sets `s_i -= 1`. The carry
`s'_{i-1} - s_{i-1}` takes 3 values (Z = 0, D = -1, W = wrap) and depends on the
cell only through `min(a_i, 2)`:

| carry \ symbol | 0 | 1 | 2 |
|---|---|---|---|
| Z | Z | W | D |
| D | W | D | D |
| W | Z | Z | Z |

Output is pointwise: Z keeps `s_i`, D gives `s_i - 1`, W gives the old `s_{i+1}`.
The entering carry must be taken from the *eventual cycle* of the composite, not
from a transient. This factorization is what makes the SSE2 core in
[`../src/snaperz_extender_fsm.h`](../src/snaperz_extender_fsm.h) possible.

**Semi-infinite bulk.** The state is always `P_t W_t W_t ...`, and one sweep at
most doubles the spatial period, since the 3-state map admits only fixed points
and 2-cycles. Doubling times: period 2 at t = 1, 4 at 2, 8 at 6, 16 at 16,
32 at 452, **64 at 246,200** — and no further doubling through 1e8 sweeps.

**Shielding is 2-adic and saturates.** `s_n` = 0 (n odd), 1 (n = 2 mod 4),
5 (n = 4 mod 8), 15 (n = 8 mod 16), **435 (16 | n)** — 0 mismatches for
n = 2..80, and independent of the tail. That uniformity is what would make a
proof possible. The **shielded-clock law**
`tau_n(y) = s_n + 2n - 1 + (A_m^{s_n}(y))_0` has 0 mismatches for n = 4, 8, 16
over all tails m <= 8. At sweep 435 a length-16q system is exactly q copies of
`P_16 = (2,0,4,0,0,1,1,0,2,0,2,0,2,1,1,0)`; at 436 a travelling carry front
starts. **The law fails at n = 32**, and since `s_n` saturates, the hierarchy is
a bounded-time phenomenon and cannot yield a fast algorithm.

**Interface fact.** One sweep of cells 0..15 depends on the whole tail *only*
through `a_16` (111 classes, one successor each), and the prefix influences the
tail only through the 3-state carry. So the two-way interface is exactly
(`a_16` in, carry out).

**112-block Toeplitz template.** For m = 15..30, `112 | q_m`, and positions
0..110 of each 112-block are universal, independent of m, with gap frequencies
`{2:56, 4:28, 6:14, 10:8, 12:4, 14:1}` summing to 450. Only position 111 (the
"hole") varies, so `p_L = 450 * (q_{L-1} / 112) + sum(holes)`.

**The renormalized map `R_r`, which supplies the states the template never did.**
Write `S_{t,r}(c,v) = (c, A_{r-1}^t v)` for the root-preserving forest clock.
Then one restored 16-cell prefix induces the tail map

    R_r = H_r . S_{30,r} . A_r^435,   i.e.   R_r(y) = J_r^c(0, A_{r-1}^31 v)
    where (c, v) = A_r^435(y),

because the reset deletes the root `(c,v) -> (0, A_{r-1} v)`, runs 30 sweeps with
the lead coordinate zero, then reinjects `c` times, for elapsed time
`435 + 1 + 30 + c = 466 + c`. Hence

    A_{16+r}^{466+c}(E_16 || y) = E_16 || R_r(y),   H_{15+r}^112(E_15 || y) = E_15 || R_r(y),

and with `y_0 = E_r`, `y_{j+1} = R_r(y_j)`, `M_r` the period,

    p_{16+r} = sum_j (466 + c_j) = 466 M_r + sum_j z(y_j),

using `z(R_r(y)) = (A_r^435 y)_0`. Verified: the closed form against the literal
112-returns for all 625 Catalan tails to r = 7; the recurrence against direct
simulation for L = 17..34; and `M_r` against the `B_L` table to L = 42.

`M_r` for r = 1..26: 1, 1, 2, 1, 4, 2, 4, 2, 8, 10, 12, 16, 24, 10, 246, 14, 78,
19, 430, 80, 824, 34, 1506, 46, 816, 395.

This *is* the 112 template: the block count is exactly `M_r` and the hole is
exactly `16 + z(y_j)`, so `450 M + sum(holes)` and `466 M + sum z` are the same
sum. The template gave weights and never states; `R_r` is the missing half. What
it does not give is a speedup — `M_r ~ p_{16+r}/469`, so the reduced orbit is
just as long, and a literal evaluation of `R_r` does the same weighted work on
`L - 16` cells instead of `L`. The ceiling is `L/(L - 16)` (2.07x at L = 31,
1.30x at L = 69, 1.016x at L = 1024), which tends to 1 and so is not a
constant-factor win. It is nonetheless a better macro than the checkpoint map,
where reducing the iteration count raised the per-iteration work by almost
exactly the inverse amount; here the 16 cells really do leave the computation.

**The obstruction to a bounded block interface.** The tempting simplification
`R_r = H_r . A_r^435` drops `S_{30,r}` and happens to hold through r = 5, then
fails at r = 6 (31 of 132 states; 151/429 at r = 7). The reason is exact: the
root-one fibre contains all phases `(1, A_{r-1}^t E_{r-1})`, at r = 6 the
canonical `A_5`-cycle has `p_5 = 13`, and `30 = 4 (mod 13)`, so the forest clock
is nontrivial. Any exact block interface recording only the root — or any
bounded alphabet independent of r — must therefore fail; the interface has to
carry an arithmetic phase. This is a rigorous reason the static-atom and
root-only cocycle approaches could not close.

**Small-tail phase map, and why phase coordinates die.** `G_m = Theta_m . H_m`
preserves root degree, so it has a root-one section
`G_{s+1}(1, v) = (1, Gamma_s(v))` with `Gamma_s = Theta_s . S_{1,s} . A_s`
(0 failures, s <= 8). This isolates the phase question inside `C_s` without
simulating a 65-71 cell extender, and the answer is bad. Representing a state as
(`A`-cycle identifier, phase mod its period) needs `Gamma_s` to carry the
canonical `A_s`-cycle into few cycles with an affine action on phase. Instead it
**scatters**: the `p_s` points of the canonical cycle land on

    s        2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17
    p_s      2    5    6   13   14   28   30   31   32   64   34  129  198  225  466  467
    targets  1    1    3    5    8   14   20   23   26   49   30  112  175  205  426  444

distinct target cycles, with target periods running up to 7,162 at s = 17. The
targets-per-point ratio climbs monotonically — 0.50, 0.67, 0.81, 0.87, 0.91,
0.95 at s = 4, 8, 10, 13, 15, 17 — so the "cycle identifier" is asymptotically as
big as the state itself, not a small label. The two cases that do land on a
single cycle settle it further: at s = 3 the induced phase map is not even
affine. This is the failure mode the phase-coordinate route was warned about, and
it happens at the smallest sizes, not asymptotically.

### The recursive lift

**The last-peak lemma is now complete.** In rightmost-peak (ECO) coordinates
`lam_n(t, r)`, `0 <= r <= h(t)`, where `h(t) = n - max{i : t_i > 0}` is the final
descent and the peak is inserted leaving `r` old down-steps after it, the final
descent takes a Motzkin step: `h(A_n t) - h(t) in {-1, 0, +1}`, decided by the
value `c = u_j` (`u = P_n(t)`, the gates strictly before the last positive
coordinate, `j = n - h(t)`) as `c = 0 -> +1`, `c = 1 -> 0`, `c >= 2 -> -1`. The
up and down counts are each exactly `C_{n-1}`. `A_{n+1}` then acts on `(t, r)`
by shifting the fibre, except for the single exceptional edge

    pi(A_{n+1}(lam_n(t, h(t)))) = P_n(t)   (not A_n(t))

on a down edge, which is the case the old lemma excluded. Zero failures over all
of `C_10`.

**`Xi` and the return recurrence.** The section
`eta_m(v) = lam_{m+1}(lam_m(v, h(v)), h(v))` is exactly the set of states whose
rightmost nonzero coordinate is 1, is followed by a zero, and has a left
neighbour `>= 2`; it has `C_m` elements. The first return `Xi_m` to it is a
permutation of `C_m` with weights `w_m`, and

    p_{m+2} = sum over the canonical Xi-cycle of w_m,

anchored at `eta_m(R_m) = A_L(R_L) = (L-1, 1, 0, ..., 0)`. The canonical weight
word is a **palindrome** for every L = 3..46 and for L = 69; that is specific to
the canonical cycle (at m = 3 only one of the three `Xi_3`-cycles is
palindromic). It is a sharply stated theorem candidate with no proof yet.

**Two-slot monodromy.** Root a base `A_n`-cycle at a state of minimum final
descent; the cyclic height word is then a Motzkin excursion. Every fibre
position `d >= 2` returns to itself after one lap, and the one-lap monodromy on
`{0, 1}` is one of the seven partial injections of a two-element set. Hence
every hole exits in under three base periods, `0 < w < 3p`, which is a property
of the ECO fibre dynamics and not of Snaperz. One scan of a base cycle therefore
yields all of `Xi` above it, and the cycle-length multiset of `A_{n+1}` follows
from the cycles of `A_n` — verified against direct enumeration through `A_15`.

**One sweep is a short recursion on plane trees.** Append an implicit zero to a
Catalan state and read it as the preorder outdegree sequence of a rooted plane
tree `T = (T_1, ..., T_d)`. Writing `U (+) V` for "append V as the last child of
the root of U", the whole sweep is

    F(.)                    = .
    F((S))                  = (., F(S_1), ..., F(S_k))
    F((T_1, ..., T_d))      = (F(T_1 (+) T_2), F(T_3), ..., F(T_d)),  d >= 2

with an equally short inverse (a leading recovered leaf can only have come from
the unary case; otherwise split the last child off the first child). This is a
structural proof that `A_n` is a permutation, and it is the exact reason
firewalls exist: after the root rewrite each child occupies a contiguous
preorder interval ending at a leaf, so the boundary gate is a no-op and the rest
of the sweep factors. The interface law is
`F(U (+) V) = F(U) (+) F(V)` once `deg U >= 2`, so the *only* cross-subtree
carry is a reversible push down a unary spine. Each local gate moves one
subtree endpoint, `ord(T_i) = lcm(1, ..., n-i)`. Zero failures over all Catalan
states to n = 11 and all pairs to combined size 10.

**Why that `lcm` is not deep, and carries no arithmetic** — worth stating,
because `lcm(1..k) = e^{psi(k)}` looks like the one arithmetic object in the
file and invites importing machinery. It has a one-line proof. Gate `i`
preserves `s = x_i + x_{i+1}` and touches only that pair. The prefix condition
at `k = i-1` gives `P = sum_{j<i} x_j >= i`, so the constraint `x_i >= i+1-P`
has lower bound at most 1; the admissible `x_i` are therefore `{1..s}` or
`{0..s}`, and `x_i -> x_i - 1` for `x_i >= 2` with `1 -> s` *is* the `s`-cycle
on `{1..s}`, with 0 fixed. Achievable pair sums are exactly `1..n-i` (and
`2..n` at `i = 0`, since `x_0 + x_1 >= 2` forces no 1-cycles there). So every
gate is a disjoint union of pair-sum rotations, one per
(prefix, suffix, `s`) class, and the `lcm` is just "an interval of rotation
lengths". Verified as stated — cycles are pair-sum rotations in every class,
n = 4..9, all `i`. Consequently `ord(A_n)` has no clean form (`~10^135` at
n = 12; n = 7 is `2^4 * 3^3 * 5 * 7^2 * 11 * 13 * 17 * 31`) and `p_n` is not
arithmetic in `n` — `p_5 = 13`, `p_9 = 31`, `p_12 = 2 * 17` all have a prime
factor above `n`, and `p_n | lcm(1..n)` fails from n = 5 on.
Measured in `research/audit/r28_gatearith.py`.

**Even lengths are 2-Dyck states.** Encoding the good pairs by
`(0, 2d) -> (0, d)` and `(2c-1, 2d+1) -> (c, d)` is a bijection from the aligned
sector onto the Fuss-Catalan states `sum b_i = n`, `2 sum_{i<=k} b_i >= k+1`,
counted by `C(3n, n)/(2n+1)`. `A` alternates the aligned and staggered sectors,
so `K_n = A_{2n}^2` is a permutation of them, and since the physical final
descent has opposite parity in the two sectors, `h(Ax) - h(x) in {-1, +1}`:
**there are no flat height steps in the even sector**. One half-sweep is a
sequential transducer with a *one-bit* carry over the pair coordinates — the
even-sector specialization of the 3-state carry automaton in the engineering
section, with the wrap case excluded by parity. The corresponding ternary ECO
lift has 8 edge types, a three-slot monodromy and `0 < w < 4p`.

**No exact merging quotient exists.** If `q(Fx) = f(q(x))` is deterministic and
some observable recoverable from `q` singles out one point of a cycle, then `q`
is injective on that cycle. On the canonical `A_L`-cycle `a_0 = L` holds only at
`R_L`, so *every* exact deterministic quotient that can still recognize
retraction has at least `p_L` states on the orbit; same for `H_m` and the
checkpoint weight `b_0 + 1`. This is a one-line argument and it retroactively
explains the whole family of failed labels below — parity, decorated parity,
defect masks, phase labels, block alphabets. What it does not rule out is an
*injective* representation with cheap random access, which is what the remaining
routes have to be.

**The pointed reversor `E_n`, and the two reflection conjectures it closes.**
There is an explicit, recursively computable involution `E_n : C_n -> C_n` with

    E_n^2 = id,     E_n A_n = A_n^{-1} E_n,     h(E_n x) = h(x).

It is built from the complete ECO lift, by induction on the fibre coordinate:

    E_{n+1}(t, r) = ( E_n(A_n t),     r    )   if r <= h(A_n t),
                    ( E_n(Phat_n(t)), h(t) )   if r = h(t), h(A_n t) = h(t) - 1

with `Phat_n(t)` the gates strictly before the last positive coordinate — the
exceptional top-fibre edge the original last-peak lemma excluded. The two cases
are exhaustive, and the fibre coordinate is never touched, so height is
preserved and `E_n` costs `O(n^2)` with no reference to `p_L` or to phase.

Height preservation is the whole content: `R_n` is the unique state of height
`n`, so `E_n(R_n) = R_n`, hence `E_n(A_n^s R_n) = A_n^{-s} R_n`. That gives, as
theorems rather than measurements:

- **global height reflection** `delta_t = -delta_{p_L - 2L + 1 - t}`, and in the
  stronger state-level form `E_L(A_L^s R_L) = A_L^{-s} R_L`;
- **the canonical `Xi`-weight palindrome** `w_i = w_{q-1-i}`, via `e_s = e_{-s}`
  on the exceptional-edge indicator — the canonical cycle self-pairs because
  `R_n` is fixed, which is also why *generic* `Xi`-cycles need not be
  palindromic;
- **dihedral pairing of every cycle** of `A_n`, not just the canonical one.

Audited independently in `research/audit/r26_reversor.py`, rebuilt from `core.A`
and the round-18 `lam`/`unlam`/`P` rather than transcribed: exhaustive over every
Catalan state to n = 12 (208,012 states), 0 failures on closure, involution,
conjugation, height preservation, `E(R) = R`, and canonical reflection. Cycle
pairing at n = 12 is 3,452 cycles, 2,702 self-reversing, 750 in 375 reversed
pairs. Compare the "global palindrome involution" dead end, which this does
*not* contradict — see that entry.

**`eps_L` is a midpoint observable.** Reflection pairs every non-fixed phase, so
the winner bit is decided at the one other fixed point of `E_L`:

    eps_L = 1  iff  p_L is even and h(A_L^{p_L/2} R_L) = 1,

equivalently `T_L = 2 p_L - L + 1` in that case and `p_L - L + 1` otherwise. An
odd canonical period can therefore never have `eps_L = 1`. Checked against the
community `T_L` table for L = 2..46: 45 of 45 lengths, 0 mismatches, and the
criterion independently regenerates `eps_L = 1` exactly at
L = 2, 4, 6, 12, 20, 21, 41.

## Data

`B_L` for L = 4..46:

    3, 9, 9, 22, 23, 23, 23, 54, 23, 117, 185, 211, 451, 451, 451, 918, 451,
    1854, 919, 1853, 919, 3729, 4667, 5608, 7493, 11251, 4673, 115288, 6559,
    36504, 8887, 201627, 37551, 386120, 15927, 705934, 21567, 382512, 185403,
    2232272, 38513, 729845, 33767

Plus `B_60 = 7,640,649`, `B_64 = 22,789,705`, `B_66 = 74,827,985`,
**`B_69 = 1,124,943`**, and `p_49 = 94,695,837`.

The gap in the even family is filled by `B_48 = 115,151`, `B_50 = 439,773`,
`B_52 = 1,527,249`, `B_54 = 5,193,987`, `B_56 = 926,585`, `B_58 = 7,763,913`
and `B_62 = 4,745,923` (`p_L = B_L + L - 1` throughout). Note `B_56 < B_54` and
`B_62 < B_58`: inside one parity the sequence is not monotone either, so a
same-parity least-squares fit is being run through a sawtooth.

**Data sources.** Everything above was computed here. There is also a community
pulse-count table maintained by the extender-hunting group, which records `T_L`
(not `B_L`) for the P = 12 column and runs much further: every `T_L` for
L <= 96, and evens through **L = 136**, 116 values in all. Export it as CSV via
`/export?format=csv` on the Google Sheets id
`1TsWICwHCjPV0JRmxDnU5vnZZ8W_uuXsN2pIdI00WaHc`.

It is worth trusting. Over the 54 lengths where it overlaps this ledger it has
**0 mismatches**, and it independently reproduces the `eps_L = 1` set
{4, 6, 12, 20, 21, 41} via `T_L = 2B_L + L - 1`. The two sources were built from
different code, so this is real mutual validation.

Two cautions. The sheet warns that its AVX loop detection can return multiples
of the true value (factors of 4, 8, 15.5, 16, 31 observed) — but that does *not*
explain the upward scatter in the P = 12 column: `L = 54` and `L = 66` sit 6.1x
and 7.4x above the even trend and both match this ledger's own C core exactly,
so that scatter is intrinsic. Second, the other eleven columns (16 gt through
56 gt period) are a different regime this project has never studied — the whole
analysis here is P = 12 — and most of their entries are `loops N`,
i.e. non-terminating.

Against a same-parity log-linear trend over 101 lengths, **L = 69 is the only
real outlier**: z = -4.99, where the next largest deviations are L = 110 at
-2.50 and L = 72 at -2.09. Its "66x below `B_66`" description understates it;
against the local trend it is ~3 orders of magnitude low and nothing else in
the extended table is comparable.

**Complete cycle censuses**, by direct enumeration of `C_n`: `A_11` has 1,218
cycles, `A_12` 3,452, `A_13` 8,978 (559 distinct lengths, max 1,874), `A_14`
24,858 (948 distinct, min 3, max 2,945), `A_15` 66,111 (1,661 distinct, min 3,
max 4,962).

**`B_L` is not smooth.** Almost every value has a prime factor above
`sqrt(B_L)`; `B_15 = 211`, `B_22 = B_24 = 919`, `B_29 = 11251`, `B_30 = 4673`,
`B_34 = 8887` are prime, and `B_69 = 3 * 374981` with 374981 prime. So `B_L` is
*not* an lcm of small independent periods, which kills the shape of every
decomposition route.

**Parity separation is robust, and the two growth slopes are now distinct.**
Odd `T_L` dwarfs even at the same `L` (`T_43 = 2,232,272` vs `T_44 = 38,513`,
58x). This entry previously said two bases were *not* supported, because
least-squares on data to L = 46 moved by ~0.15 with the window — odd 1.3575
(L >= 16) / 1.3553 (L >= 24) against even 1.2004 / 1.1762, while `B_60/64/66`
pushed the even base to 1.3471 over L = 40..66, which made a common base with
different constants look admissible.

The community `T_L` table (see Data sources) settles it. With evens to L = 136
and odds to L = 95:

| window   | parity | base       | 95% CI           | R^2  |  n |
|----------|--------|------------|------------------|------|----|
| 16..136  | even   | **1.2281** | [1.2033, 1.2534] | 0.96 | 61 |
| 16..96   | odd    | **1.3774** | [1.3442, 1.4114] | 0.95 | 40 |

Slope difference `0.1147 +- 0.0128`, **t = 9.0**; a common base is rejected. The
even base is stable in every window from L = 40 up (1.2268 +- 0.0118 over
L = 40..136) — the old 1.3471 was a short-window artifact, exactly as this entry
suspected. Fits are on `T_L`, so entries with `eps_L = 1` contribute
`2B_L + L - 1`; that is a `<= log 2` additive wobble, biasing the slope by
`<= 0.007` over these windows. **Different finite-range slopes: settled.
Different asymptotic bases: still conjectural.** Mean checkpoint gap is ~4.1 for
both parities, so odd extenders are not slow per checkpoint — they visit more
checkpoints.

The normalized exponents are suggestively close but do **not** agree: against
each sector's own state space (odd Catalan `4^L`, even Fuss-Catalan `2.598^L`),
`log(1.3774)/log 4 = 0.2310` and `log(1.2281)/log 2.598 = 0.2152`, 7% apart. If
those were equal the parity gap would be a pure state-space-size effect. They
are not, and this is the shape of numerology that has been buried twice here, so
it is a target to test, not a result.

**`E_L` is not on a distinguished orbit.** At L = 13: 742,900 states in 8,978
cycles, 559 distinct lengths, max 1874, mean 82.7, and `p_13 = 129` is
unremarkable — 40.6% of states lie on shorter cycles. Max cycle is ~`2^0.9L`
against a `2^2L` state space.

**Leftward influence** is <= 1 cell/sweep and bursty (long stall, short burst);
rightward is effectively instant. Front arrival does **not** explain `B_L`: at
L = 43 the root first differs from the bulk at 2,323 sweeps against
`B_43 = 2,232,272`.

## The L = 69 outlier

`B_69` is ~66x *shorter* than `B_66` and far below its same-parity neighbours.
No mod-16/32 family reproduces it, the "64 + 5 resonance" story is unsupported,
and L = 133 is untestable (a run passed 5e9 sweeps). It is now diagnosed
dynamically, though not derived algebraically.

Let `ell(x)` be the length of the longest good-paired even prefix. At all 2,394
16-clock returns `ell >= 56`, so violations sit in a 13-cell tail; 1,020 of them
have `ell = 68`. Neighbours L = 65..83 all lose it at outer step 781-787;
69 does not (`ell(x_780) = 56`, `ell(x_781) = 68`) and they then decay to
`ell = 42`.

Off-section, over all 268,128 checkpoints, `min ell = 0`. The only 16 sub-56
checkpoints are indices 268,111..268,126 (outer step 2,393, offsets 95-110),
collapsing monotonically 54 -> 0 into the return to `E_69`. **The defect escapes
exactly once, and that escape is the closure.**

**56+13 rigidity.** Of the 2,394 16-clock sections, 1,603 have a clean Dyck cut
at 56, and every one of those 13-cell tails lies on the 129-state canonical
`A_13` cycle (0 off-cycle). The same holds at `56 + s` for s = 9, 11, 13, 15
with `p_s` = 31, 64, 129, 225 and lock offsets `c_s` = 14, 28, 30, 30
(`phi = T + c_s`, 0 failures to j = 780). **s = 17 fails immediately** (all 780
clean tails off-cycle) — the 16-barrier again.

**Phase-reflection law.** On clean sections `z(P_-j) = z(P_j)` and
`phi_j + phi_-j = z(P_j) + 43 (mod 129)`, over 1,588/1,588 nontrivial pairs.
This is the project's first genuine partial fast-forward, but it moves the tail
*phase* only; there is no cheap `P_780 -> P_1614` prefix map, so it is not a
whole-state macro.

It holds far more widely than the 16-clock sections. Over **every** checkpoint of
the L = 69 orbit there are 179,665 clean 56-cuts, of which only 131 have their
13-tail off the canonical `A_13` cycle, and the law holds on all 178,354
reflected pairs with the same constant 43.

**Where the constant comes from, and where it does not.** `43 = 56 - 13 = n - s`.
Reading it instead as `n + B_s - 1` is an artifact: since `p_s = B_s + s - 1`,
`n + B_s - 1 = n - s (mod p_s)` identically, so `B_13 = 117` cancels and carries
no information. The law `phi_j + phi_-j = z(P_j) + n - s (mod p_s)` is **not**
general over clean `n|s` cuts. Measured over every clean cut of every orbit to
L = 44 plus L = 69:

- it holds at `n = 8` (L = 10..20, 26..32), at `n = 16` (L = 18..36, every
  s = 2..20 bar three), at `n = 24` (L = 35), `n = 32` (L = 42), and at both
  `n = 56` and `n = 48` for L = 69 — the L = 69 `48|21` cut gives 186 pairs at a
  fixed constant 27 = 48 - 21;
- it fails at generic `n`. Inside L = 69 itself the `60|9` cut has no fixed
  constant at all over the full orbit, and on the round-14 112-sections it has a
  fixed constant 25 against a predicted 20;
- and it is not universal even at the cuts that usually work: at `(n,s)` =
  (16,6), (16,8), (16,14), (8,14), (8,16) the constant is not constant, and at
  (24,4), (24,6) it is fixed but wrong.

Every cut that carries the law has `8 | n`, but that is not sufficient, and no
sharper characterization survived. Treat `43` as a measured fact about the
`56|13` cut of L = 69, not as an instance of a general law.

**The boundary-event cocycle explains it.** Cut at `n` and let the last prefix
cell drive the tail: its value at the moment its gate fires is 0, 1 or `>= 2`,
which applies `A`, `A.Q` (delete the tail root) or `A.I` (inject one unit) to
the tail. In the rotating frame `rho_t = A^{-t} T_t` the ordinary `A` sweeps
vanish entirely and only deletion and injection events move the state, so for a
tail on the canonical cycle the rotating phase `d_t = phi_t - t (mod p_s)` is
piecewise constant. At L = 69, `n = 56` the 179,534 on-cycle clean checkpoints
fall into exactly **84 constant-`d` runs separated by 83 dirty excursions**, and
erasing the no-op sweeps leaves only **four distinct event words**: length 64
(x2, jump 30), 192 (x40, jump 39), 1,487 (x40, jump 14) and 89,214 (x1, jump
123). The jump word is the palindrome

    30, (14, 39)^20, 123, (39, 14)^20, 30

with total 2,303 = 110 (mod 129). A palindromic jump word forces
`d_i + d_{m-i} = d_0 + d_m` by partial-sum cancellation, which is the defect
reflection; combined with the checkpoint-time identity
`t_k + t_{N-k} = p_L - L + z(x_k)` and `p_69 - 69 = 62 (mod 129)` it gives

    43 = 62 (checkpoint-time centre) + 110 (cocycle holonomy)   (mod 129),

the first derivation of the constant with dynamical content in it. Both
`n + B_s - 1` and `n - s` remain numerology; this is not. What is *not* proved
is why the prefix generates a palindromic 83-excursion word in the first place —
that is still an exhaustive computation.

The short and long excursions are literal subset reset words: applied to all
129 phases, the short word sends `alpha + {0,1,3,5}` to `alpha + 42` and the
long one sends `alpha + {0,1,3,5,7,22,52}` to `alpha + 19`, in both cases
because the surviving states come to differ only in their root and one deletion
event erases that. The real orbit enters the long reset at `alpha + 5`, and
after the outer step 780 that admitted set is `{18,19,21,23,25,40,70}` — which
contains the uniquely closing phase 23 and the near-misses 19, 21, 25. So the
reset word explains the *basin*; the prefix dynamics still does the final
selection.

The same 84 intervals show up as state synchronization: over 753,124 of the
1,125,011 sweeps (66.9%) cell 55 is zero and cells 56..68 sit on the canonical
`A_13` cycle, in exactly 84 maximal runs. Truncating each run to whole `11 x 129`
clocks gives the 514 complete `A_13` clocks (729,366 sweeps, 64.8%) that a
cross-dimensional parse of the L = 37 itinerary sees as one repeated tile. The
tree recursion above says why the tail is autonomous at all: the boundary gate
sits at a leaf, so the tail is a separate recursive `F`-call.

`Xi`-sections make the same structure visible at a third resolution. At L = 69
the canonical `Xi`-weight word has period 29,265 and the exact distribution

    13^64, 21^11680, 28^11680, 31^5758, 77^2, 379^40, 2597^40, 254167^1,

summing to `p_69`; the 40 + 40 + 1 rare weights are the envelopes of the 40 + 40
+ 1 dirty excursions. At this length, and only at this length among those
tested, membership of the `Xi`-section is a five-symbol local pattern: `UDDDD`
in the height word, with 0 false positives and 0 false negatives over the whole
orbit. At L = 31, 35, 37, 43, 45 the pattern still implies section membership
but catches only 5%, 13%, 1.4%, 8% and 15% of it.

**The obvious proof of the law is dead.** The natural route was to combine the
already-proven telescoping identity `t_k + t_{N-k} = p_L - L + z(x_k)` with phase
lock `phi_k = t_k + c`, which would force the constant to be `p_L - L + 2c`. But
phase lock is *not* global: over all clean checkpoints at L = 69, `n = 56` there
are 62 distinct lock offsets with the mode covering only 48.6% (round 14 saw a
single offset because it only looked at 112-sections with `j <= 780`), and even
at `n = 16` there are exactly 2 offsets at 50/50 for every L. The forced constant
would be 93, not 43. So the reflection law survives precisely where phase lock
fails, and is a genuinely separate fact rather than a corollary. It remains
unproven.

**The reflection map is not a Catalan involution.** The natural candidate is
`Psi_m = H_m . Theta_m`, which exactly satisfies `Psi_m(E_m) = E_m` and
`z(Psi_m(x)) = z(x)`, and reflects the canonical `H`-cycle through m = 5. It
fails at m = 6: with `b = H_6(E_6) = (3,0,0,1,2,0)`, `H_6 Theta_6(b) =
(5,0,0,0,0,1)` while `H_6^{-1}(E_6) = R_6 = (6,0,0,0,0,0)`. The discrepancy
again lives on the period-13 `A_5`-cycle. So `P -> P^dagger` needs a forest-phase
correction; it is not a relabelling of plane trees.

**129-phase uniqueness** is the strongest L = 69 statement to date. Holding the
real 56-cell prefix at outer step 780 and substituting all 129 tail phases:
{19, 21, 23, 25} survive 834 outer returns, {19, 21, 23} survive 1,614, and
**only phase 23 reaches `E_69`**, at exactly return 1,614. The real phase is 23.
Also confirmed: the first 56 coordinates at step 780 are identical for all odd
L = 65..83, all at `T = 366,482` sweeps; the injected root `q = (A_s^293 y)_0` is
4, 2, 6, 9 for L = 65, 67, 69, 71, so L = 71's violent escape is just `q` odd.

**Checkpoint-time reflection identity.** `t_k + t_{N-k} = p_L - L + z(x_k)`,
where `t_k` is accumulated `A`-time and `z` is taken on the full checkpoint state
with `z(E_L) = L`. 0 failures over all 268,129 indices at L = 69. The proof is a
one-line telescoping, conditional only on the gap palindrome.

**Universal terminal cone.** Final checkpoints factor as `U_r || Z_{L-e_r}` with
`Z_s = (1, s-1, 0, ...)`, `U_r` independent of `L`, and `e_{r+1} - e_r` the
112-block gap word. Verified for L = 25..43 odd and 69; every deviation is
exactly one step past the right boundary, as the law itself predicts. L = 69
first escapes at `e_16 = 54`, i.e. `U || Z_15`, so escape-equals-closure is
forced. It is `O(L)` instead of `O(L^2)`, but removes only the last `O(L)`
checkpoints.

## Rounds 27-29: what carries the L = 69 excursion palindrome

The last L = 69 conjecture is the 84-excursion palindrome that `43 = 62 + 110`
rests on. The reversor proves the canonical height and exceptional-ECO
palindromes but was not shown to reach the `56|13` defect classification. This
round locates the mechanism. Harness `research/audit/r27_reflect.c`,
`r27b_eco.c`, `r27c_bridge.c`.

Rooted at `E_L`, the reversor acts on time as the involution

    sigma(t) = (p_L - 2L + 2) - t   (mod p_L),

the same centre the height reflection uses. Three findings, in order.

**The `56|13` section is not `sigma`-invariant — not even slightly.** Of 268,128
checkpoints, **0** have their `sigma`-partner a checkpoint; of 179,534 clean
on-cycle checkpoints, **0**. So no proof can go through "the section is
preserved", which is exactly why the reversor did not close this. Nor is it an
interval reflection: excursion `e` maps onto its partner only up to an offset
that is always **odd** and bounded by 29, and 76 of 83 partner lengths differ
(by an even amount, up to 26).

**The exceptional-ECO symmetry is exact here, with the state centre.** At
L = 69 there are **29,265** exceptional edges — exactly the `Xi` period — and
`e_s = e_{C-s}` holds on **29,265 of 29,265**, while the edge-centred variant
`e_s = e_{C-1-s}` holds on **0**. Take `C = p_L - 2L + 2`. Deriving the induced
map on gaps from that invariance rather than assuming a rooting gives `pi` well
defined on all 29,265 gaps with `w_{pi(i)} = w_i` on **29,265 of 29,265** — the
`Xi`-weight palindrome, independently confirmed at L = 69 in the correct
rooting. (Assuming the `E_L` rooting instead gives a spurious 77/83; a rotated
palindrome is not a palindrome, and that artifact is easy to walk into.)

**The excursions are the rare gaps.** The multiplicities coincide exactly:

    excursion jumps :  30^2    39^40    14^40    123^1      = 83
    rare Xi weights :  77^2    379^40   2597^40  254167^1   = 83
    common weights  :  13^64 + 21^11680 + 28^11680 + 31^5758 = 29,182 = 29,265 - 83

Excursions contain almost no exceptional edges: they sit *inside* the 83 rare
gaps. The correspondence (anchored at `b0`) is an order-preserving injection,
the enclosing weight determines the jump single-valuedly
(`13/77 -> 30`, `2597 -> 14`, `379 -> 39`, `254167 -> 123`), the rare gaps are
closed under `pi` acting as `r -> 80 - r` on cyclic rank on all 83, and the
excursion pairing `e <-> 82 - e` agrees with it on **79 of 83**. Anchoring by
overlap instead, **82 of 83** excursions contain exactly one rare gap and one
contains none.

So the mechanism is identified: **the excursion palindrome is the proved
`Xi`-weight palindrome transported through an excursion/rare-gap
correspondence.** It is not yet a proof — the correspondence is exact on 79-82
of 83, and closing the residue is the remaining work. Note also this is one
length; nothing here has been checked at another `n|s` cut.

**Round 29 correction: that conclusion was overdrawn.** The rare-gap
correspondence at L = 69 is real as a correspondence, but it is *not* what
carries the palindrome. Checked at the other `56|s` cuts with round 24's run
counts (harness `research/audit/r29_cuts.c`, memory-lean rewrite that stores
event lists only, so a 73M-sweep orbit fits in 3 GB):

| L  | cut   | excursions | jump palindrome | R24 constants | rare-gap overlap |
|----|-------|-----------:|-----------------|---------------|------------------|
| 64 | 56\|8  | 933       | **fails** 383/933 | {0,6,18,20} | 75/933 unique    |
| 68 | 56\|12 | 47        | 47/47           | {10}          | **0/47** (46 none) |
| 69 | 56\|13 | 83        | 83/83           | {43}          | 82/83            |
| 70 | 56\|14 | 229       | 229/229         | {42}          | 0/229            |

At L = 68 a threshold even exists making rare-gap and excursion *counts* match
(47 = 47 at w > 11207) and the overlap test still returns zero — a count match
alone is worthless, which is Rule 6's cousin. The L = 69 containment is an
accident of its extreme gap structure.

**The actual mechanism is a defect-sum law, and it is a two-line derivation.**
Subtract the round-14 time identity `t_k + t_{N-k} = p_L - L + z(x_k)` from the
round-15 reflection law `phi_k + phi_{N-k} = z(x_k) + C` (single constant):

    d_k + d_{N-k}  =  C - (p_L - L)   (mod p_s),   d_k = phi_k - t_k,

constant over clean reflected pairs, k the *global* checkpoint index. A
constant defect-sum maps constant-`d` runs to constant-`d` runs reversed,
boundaries to boundaries, and reflected jumps satisfy
`(K - d_e) - (K - d_{e+1}) = jm[e]` — the jump palindrome. Measured:

| L  | predicted `C - (p_L - L)` | observed mode        | off-mode pairs |
|----|---------------------------|----------------------|----------------|
| 68 | `10 - 14 ≡ 30` (mod 34)   | 30 on 176,586/176,587 | the `k=0` self-pair |
| 69 | `43 - 62 ≡ 110` (mod 129) | 110 on 178,354/178,355 | the `k=0` self-pair |
| 70 | `42 - 74 ≡ 166` (mod 198) | 166 on 185,166/185,167 | the `k=0` self-pair |
| 64 | no single `C` exists      | 5 values, mode 24 at 99.992% | 77 real violations |

The L = 69 sum 110 is exactly the round-16 holonomy (`43 = 62 + 110` was this
law read backwards). And L = 64 is the exact contrapositive: 77 violated pairs
out of 987K scramble 550 of 933 palindrome positions, because one violation
inside a run shifts every later boundary — the palindrome is brittle exactly as
the derivation predicts.

**Status change.** The 84-excursion palindrome at L = 69 is no longer an
independent conjecture: it is a corollary of (a) the round-15 single-constant
phase-reflection law (still open) and (b) the round-14 time identity (proved by
telescoping, conditional on the checkpoint-gap palindrome). Loose ends kept
honest: excursion counts match round 24 exactly at 68/69/70 but give 934 runs
vs their 937 at L = 64 (definition wrinkle at off-cycle interruptions, not
chased). Also confirmed in passing: the exceptional-edge set is
`sigma`-invariant at every length tested, now including 12,119,831 of
12,119,831 edges at L = 70 — a large-scale confirmation of the reversor
theorem.

**The full `56|s` census, and the two failure modes.** Every even length with a
walkable orbit, same harness:

| L  | s  | `p_s` | on-cycle clean | excursions | defect-sum values | palindrome |
|----|----|------:|---------------:|-----------:|-------------------|------------|
| 58 | 2  | 2     | 364,853        | 0          | —                 | trivial    |
| 60 | 4  | 6     | 115            | 25         | 1 (`0`)           | 25/25      |
| 62 | 6  | 14    | 32,221         | 356        | 1 (`2`)           | **346/356** |
| 64 | 8  | 30    | 993,409        | 933        | **5**             | **383/933** |
| 66 | 10 | 32    | 187,557        | 560        | **4**             | **230/560** |
| 68 | 12 | 34    | 177,224        | 47         | 1 (`30`)          | 47/47      |
| 69 | 13 | 129   | 179,534        | 83         | 1 (`110`)         | 83/83      |
| 70 | 14 | 198   | 188,456        | 229        | 1 (`166`)         | 229/229    |
| 72 | 16 | 466   | 177,976        | 57         | 1 (`266`)         | 57/57      |
| 74 | 18 | 468   | 8              | 1          | —                 | trivial    |

Two distinct failure modes, both now observed:

- **Constant multiplicity** (L = 64, 66): the defect-sum takes 4-5 values and
  the palindrome shatters. The failing tails are exactly s = 8, 10
  (`p_s` = 30, 32); s = 4, 6, 12, 13, 14, 16 all carry a single constant.
- **Orphan closure** (L = 62): a single constant is *not* sufficient. There the
  off-cycle clean count dwarfs the on-cycle (1,084,430 vs 32,221, ~24% of
  clean checkpoints orphaned), run boundaries land on orphaned positions, and
  the palindrome fails 346/356 despite a perfectly single-valued defect-sum on
  the pairs that do exist. The derivation therefore needs *both* hypotheses:
  a unique constant, and reflection-closure of the clean set at run boundaries.

Also new here: **L = 72 at `56|16` is completely regular** (single constant 266
mod 466, palindrome 57/57) — round 24 saw L = 72 as a low-locking mess, but it
was looking through the `64|8` cut; the cut was messy, not the length. And the
s = 16 tail behaves despite the 16-barrier.

The open problem, sharpened: prove the single-constant reflection law where it
holds, and explain why the tails s = 8, 10 (`p_s` = 30, 32) admit several
constants while s = 4, 6, 12, 13, 14, 16 admit one.

## Dead ends

Do not retry these without genuinely new information.

- **Naive decreasing potential** — invalid for the real dynamics.
- **Explicit marker coordinate** — correct, but the same `O(L * B_L)` and 41x
  slower than the C core. Its only merit is asserting `r in {0, L}` at retraction.
- **Simple odometer** — loses information, and the carry rate is a flat ~12% for
  all L in 8..28 with no digit hierarchy, so the recursion costs *more* than
  direct simulation.
- **Fitted numerical formulas** — counterexamples.
- **Parity quotient** — not Markovian. The **decorated** parity quotient is
  exact, and in fact an invariant subsystem rather than merely a quotient: good
  pairs `(0, 2q)` / `(2p+1, 2q+1)`, `A_{2n}` alternates the aligned sector `S0`
  and staggered sector `S1`, `S0` is in bijection with full ternary trees
  (`|S0_n| = C(3n,n)/(2n+1)`), even checkpoints are ordered *pairs* of ternary
  trees with count `N_n = 2/(3n-1) * C(3n-1,n-1)` = 1, 2, 7, 30, 143, 728, 3876,
  21318, 120175, and `tau = b_0 + 1 = 2(1 + ldepth(P))`. It is still `2.598^L`.
  Worse, **state-space counts measure the wrong thing**: the orbit is
  exponentially shorter than the space (`p_49 ~ 9.5e7` against `2e20`).
- **Nested first-return recursion `F_k`** — telescopes to the flat orbit, zero
  cache hits.
- **Suffix-independent stems** — exactly k = 1, 3, 7, 8, 16 with times 1, 2, 8,
  8, 112, swept to 45. The law `H^112(1^16 u) = 1^16 J(u)` is real, and `J` does
  have a closed form: it is `R_r` above. But "closed form" here means an explicit
  466-sweep composition, not a cheap one, so Rule 1 still bites.
- **Recursive 16-cell peeling** — the one route the `R_r` derivation suggested,
  and it dies at the first step. If the renormalized tail dynamics re-entered the
  shielded form we could peel 16, then 32, then 48 cells; instead `R_r(E_r)` never
  has an `E_16` prefix for any r = 17..30, and over the full `R_r` orbit the
  shielded set is essentially never revisited (0 visits for r = 17..24 and 26; 4
  of 816 at r = 25). The longest run of leading ones anywhere on those orbits is
  1-7 cells. So the 16-peel is a one-off, the `L/(L-16)` ceiling stands, and
  there is no hierarchy to recurse on.
- **Proving the phase-reflection law from phase lock** — see the L = 69 section.
  Phase lock is not global, and the constant it would force is wrong.
- **Phase coordinates, i.e. representing a state as (`A`-cycle identifier, phase
  mod its period)** — `Gamma_s` scatters one source cycle across ~`p_s` target
  cycles already at s = 17, and is not affine even when it does not scatter.
- **Clean 4 -> 8 -> 16 -> 32 shielding hierarchy** — fails at 32. Related dead
  guesses: `s_n = B_n - n` (predicts `s_32 = 6527`, measured 435) and growing
  shielding (it saturates).
- **Bounded-width front / constant propagation speed / front arrival explains
  `B_L`** — all false.
- **Box-ball systems and solitons** — rightward propagation is unbounded (a
  3-cell perturbation reaches the end of a 50,100-cell array at t = 2), the far
  tail flips to a *different* vacuum word `(0,1,2)` rather than a shifted
  `(1,0,2)`, and defects never re-localize — they merge instead of passing
  through. Mechanism: the mean-one vacuum never resets the carrier.
- **Yang-Baxter / transvections** — real algebra, no explanatory value.
- **`Xi`-weight multiplicity counting** — the idea was that `p_L = sum w` makes
  the period a dot product of a *weight alphabet* against a multiplicity vector,
  which is a counting problem rather than a simulation problem and so is not
  touched by the no-merging theorem. It dies on the alphabet size. Distinct
  `Xi`-weights `d(L)` over L = 16..46 fit `d ~ p_L^0.486` (R^2 = 0.80) against
  R^2 = 0.46 for `1.107^L` and 0.48 for `L^3.10` — real compression, still
  exponential, and the same signature as the LZ78 grammar's `p^0.68`. Rule 6
  again. The seductive datum was L = 69's 8-symbol alphabet
  (`13^64 21^11680 28^11680 31^5758 77^2 379^40 2597^40 254167^1`, summing
  exactly to `p_69 = 1,125,011`); L = 43 has 271 distinct weights at a
  comparable period. **That was generalizing from the one length this ledger
  documents as unrepresentative** — a mistake worth naming, because the L = 69
  distributions are the prettiest objects in the file and will tempt again.
  Related: `Xi_period ~ p_L^0.818`, so the section itself only compresses the
  orbit by `p^0.18`. Swept via `research/build/r18_xi L --dist`, L = 4..46.
- **Global palindrome involution `J A J = A^{-1}`** — vacuous *as stated*, since
  every permutation is strongly real; and preserving every prefix-reset set needs
  one common centre, which fails (at L = 16, `S_L` forces c = 0 but `S_1` gives
  450, and `S_9` has no centre at all). The palindromes are real but live in
  *induced* time, one involution per level: per-level centres at L = 16,
  k = 1..8 are 450, 448, 448, 442, 428, 428, 428, 428, with gaps `2^j - 2`.
  **Do not use this entry to bury the reversor `E_n` (Settled theory).** That
  map is not this claim: the content there is *height preservation*, which is
  what forces `E_n(R_n) = R_n` and makes the reflection non-trivial. Bare
  existence of a conjugating involution is what is vacuous; a height-preserving
  one is not, and it is proved.
- **lcm / factorization / OEIS lookups.**
- **Sparse perturbation of rowmotion** — rotating-frame states `O^{-t} A^t(E)`
  are essentially all distinct (115,314 of 115,318 at L = 31).
- **Polynomial-dimensional `F_2` linearization** — the linear complexity of
  `(A^t E)_0 mod 2` is 3,753 / 11,278 / 115,318 / 201,661 for L = 25/29/31/35,
  i.e. exactly `p_L`.
- **Affine-region skipping** — branch vectors are 86-95% distinct with **0**
  consecutive repeats.
- **Independent recursion over root subtrees** — `G_m` fails to preserve the
  root-subtree-size multiset on 8,290 of 16,796 states at m = 10.
- **Small fixed-label permutation lift.**
- **Finite look-and-say atoms** — 52,867 distinct 16-cell block states by sweep
  12,000 at q = 32, still growing.
- **Cyclotomic node budget / sandpile SNF fingerprint** — numerology;
  `budget <= m` fails at m = 36, 38, 41, 42, and large primes in `q_m` have
  `ord_2 = Theta(p)`.
- **Tropical / min-plus repeated squaring** — the sweep is not monotone:
  `s = (1,6) -> (6,6)` but `(2,6) -> (1,6)`.
- **Transducer composition** — the anti-diagonal count `N(t)` grows at least
  `t^4` with the exponent *rising* with L, so `N(1e7) ~ 1e28`; and
  `N(t) <= B_L * L` means the transducer is the size of the orbit it replaces.
  Closed by scale, not by growth law.
- **Bulk-period <-> length correspondence** — stage 32 spans t in [452, 246200]
  and contains closures for 20 lengths with `B_L` in 918..201,627; L = 19 closes
  at period 1.68x its length, L = 40 at 0.80x.
- **NP / ETH hardness** — the canonical input is unary `L -> T_L`, so it is
  sparse/tally and hardness would collapse classes. A black-box `Omega(T_L)`
  bound holds, but the rule is explicit, so it proves nothing.
- **Compressing the L = 69 separator word** — linear complexity 2,390-2,394 out
  of 2,394.
- **Polynomial-size grammars for the itinerary** — the exact LZ78 phrase count
  `g` of the final-descent `U/F/D` word is a real compression (0.47% of the word
  at L = 66, well under the `1/log N` a random word would give) but it is *not*
  polynomial in L. Over the whole even family L = 30..66, `g/L^3` ranges from
  0.0099 at L = 44 to 1.2235 at L = 66 — a 124x spread, and not monotone.
  Regressing gives `log g = 0.684 log p + 0.17` with `R^2 = 0.956`, against
  `R^2 = 0.80` for the best power of L; that is `g ~ p^0.68 ~ 1.19^L`.
  Exponential with a smaller base is still exponential. The `g ~ 1.224 L^3`
  reading is a one-point fit at the largest measured L, which is also the
  largest `p` — rule 3 in a new costume.
- **2-automatic coordinate formulas** — the retraction indicator has exactly one
  1 per period, so its 2-kernel is at least the odd part of `p_L` (57,659 at
  L = 31, and all of `p_L` when `p_L` is odd). The failure of `F_2` linear
  recurrences was not the linearity.
- **Finite-state wreath recursion on plane trees** — the natural one-hole
  version is dead by measurement: with branch gadgets `F` and `C` and two probe
  forests, all `2^d` contexts of depth `d` give distinct output-degree responses
  for every `d <= 15`, and the residual-section count over one-hole contexts of
  size `h` grows 1, 2, 4, 8, 17, 39, 102, 295, 912, 2,945, 9,798 (these two
  censuses are reported, not independently reproduced here). What survives is
  only the narrower hope that the canonical orbit visits a small *unary-spine*
  sub-language of contexts.
- **Odd = even ternary system plus O(1) defects** — with `delta(x)` the minimum
  number of bad pairs over all choices of one unmatched even cell, the canonical
  orbit reaches `max delta = (L-15)/2` exactly, for every odd L = 17..45. The
  defect count is extensive: it is a defect gas, not one moving particle. (The
  first state attaining the maximum always has an exact 16-cell staggered even
  core and a wholly defective suffix — the 16-barrier from a new direction. At
  L = 69 `max delta = 7`, mean 1.4478, and 23.8% of sweeps have `delta = 0`,
  which is the 56-cell rigidity again.)
- **Low-degree polynomial invariants** — through degree 3 the invariant space is
  reported to be exactly `1, S, S^2, S^3` in the total mass `S`. Finite evidence
  only, and not reproduced here.
- **Expanding the recursive lift on demand** — following `Xi`/`Omega` and fully
  expanding every touched lower-dimensional cycle is reported to cost 1.5x,
  4.5x, 12.9x, 24.9x, 67.5x the target period at physical L = 22, 24, 28, 30,
  32. The lift is cheap only in batch over a base cycle, and the canonical orbit
  does not stay over few base cycles.

An unverified but plausible reframing, with no algorithm attached: read the
sweep as a **"frozen-top row whirl"** on staircase partitions (order ideals of
the type-A root poset), with local map `(L, ..., U-1)(U)` — rowmotion with each
row's cycle top cut out. This explains why the rowmotion and promotion
comparisons failed. Made precise, it is the lazy-winching form
`y_j = L - p_{L-1-j} + j` swept `j = L-1..1`, which is exact. Note that
`L_j = W_j . S_j` is **false** as usually stated — 15.79% of gate applications at
L = 10 violate it, in two equal halves (at the ceiling `y_j = 2j`, and where the
ceiling binds but `S` fires anyway). It becomes exactly true with the guard
*apply the swap only when the staircase ceiling is slack*, `2j >= y_{j+1} - 1`
(0 failures over all of `C_L`, L = 5..11). The useful part is the defect
*density* — 31-38% of gates, 7.6-13.3 per sweep — which is what kills the
sparse-perturbation hope.

## Rules for judging new proposals

1. **A return or macro map only speeds anything up if its successor state can be
   produced in `o(weight)` work.** Check that before benchmarking. The same
   obstruction has been rediscovered six times: the 112 template, `G_16`, `rho`,
   the "4.2x checkpoint speedup", the "L=18 gap word", and the 16-clock
   renormalized map `R_r`. The one thing that did survive the rule is `R_r`'s
   explicit form, and only because it drops 16 cells rather than reducing the
   iteration count — worth `L/(L-16)`, which tends to 1. Measured
   flat-vs-checkpoint gate counts at
   L = 12..24 are 1.19, 1.09, 1.08, 1.07, 1.06, 1.06 — converging to 1. **Any
   summary calling the checkpoint map a constant-factor win is wrong.**
2. **Repeated squaring needs `|A^t|` polylog in t, not polynomial.** An earlier
   "polynomial implies squaring wins" framing was simply wrong.
3. **Sampled orbit statistics never converge here.** 4x fewer sweeps halves
   `N(t)` at every t, because a longer window contains genuinely more patterns.
   This artifact produced a false subexponential reading four separate times.
4. **Before believing any hierarchical or odometer proposal, measure whether the
   carry rate decays with depth.** It does not — it is flat at ~12%.
5. **External research ledgers reliably drop the closed routes and the
   deliverable.** When checking one, confirm that tropical, the `N(t)` scale
   argument, the bulk-period correspondence, the 112 template, the ~2% odds, and
   the AVX2 work are all present.
6. **Never accept a `poly(L)` reading of a quantity that was only measured at a
   handful of lengths — regress it against `p_L` as well.** Everything here that
   looks polynomial in L at the largest measured length has turned out to be a
   power of `p_L`. The check costs one extra fit and has now caught the same
   error twice (`N(t)`, and the LZ78 grammar size).
7. **A compressed representation that has to be built by reading the orbit has
   not beaten anything.** Answering a query over a grammar in `o(w)` symbol
   operations is only a speedup if the grammar itself can be constructed in
   `o(p_L)`. This is rule 1 wearing a different hat, and it is easy to miss
   because the measured query speedup is real.

## Engineering, which is the actual deliverable

- The 3-state automaton is a real constant-factor win, but the **cumulative
  rewrite alone is a pessimisation** — the win is entirely the factorization it
  enables. At L = 49 (`p = 94,695,837`): existing core 35.1s, cumulative
  branchless scalar 27.9s, cumulative + FSM + SSE **7.25s**.
- Against this repo's linked-list baseline, which already skips empty segments,
  the real win is only **~1.5x and on the fallback path alone**. It does not help
  the AVX2 path, structurally: the wavefront parallelizes *across* pulses, so
  there is no within-pulse serial carry to collapse, and it already keeps the
  prefix sum. L = 49: fallback 21.7-23.3s, FSM 14.8-15.7s, AVX2 1.4-1.7s.
- **Vector width is irrelevant to the FSM**; only automaton chain length matters
  (SSE2 15.69s vs AVX2 15.64s; 8 symbols per lookup with a 1.5 MB table gives
  11.18s). The 4-symbol / 3 KB path was kept, since that path targets
  small-cache pre-AVX2 CPUs.
- Keeping the simulation state in registers is a measured **2.1x** on AVX2
  (L = 49: 2.20s -> 1.03s).
- `src/blockless.h` implements the good-pair reduction and matches
  `research/snaperz_p12.py` exactly (L = 49 gives `p = 94,695,837`; L = 31 gives
  `B`, L = 41 gives `2B + L - 1`). Good cross-validation of both.
- **The blockless wavefront's cost does not depend on L.** Measured 9.5-12 ns
  per sweep, flat over L = 21..64. It holds the whole extender in lanes, so the
  only thing that matters is whether the extender *fits*. This is also the
  reason `R_r` is useless to it: `R_r` removes 16 cells and leaves the sweep
  count alone, and the wavefront does not care about cells.
- **What actually cost time was the fit, not the length.** Two 32-lane uint8
  windows hold `kSegCount <= 64`, so L >= 65 silently fell back to the scalar
  sweep — about 30x slower per sweep. `src/blockless_avx512.h` widens the
  windows to 64 lanes and carries L <= 127; the ceiling is 127 rather than 128
  because lanes are compared as signed bytes and `R_L = (L, 0, ...)` has to fit
  in one. **L = 66: 24.9s -> 1.29s (19x). L = 69: 0.41s -> 0.018s (23x).**
  Validated by `research/audit/r15h_avx512.cpp`: 30 lengths cross-checked
  against the 256-bit windows with 0 mismatches, plus every tabled `B_L`.
  The wide path is *slower* where both apply (0.7x at L = 64 — the byte shift
  costs a cross-lane permute where 256-bit gets one `alignr`, and the wider code
  clocks lower), so it is tried second and only where the narrow one does not
  fit. Against 1.5x from `R_r` on the same lengths, the fit is worth ~13x more
  than the renormalization, and needs none of the theory.
- The `R_r` renormalization does give a real win on the *scalar* path, where
  cost is `O(p_L * L)`: measured 1.2-3.1x over L = 31..69, tracking the
  predicted `(L-1)/(L-17)` and often beating it on cache. It is also a correct
  counter rather than just a period formula — `eps_L`, and so `T_L`, come out
  right, which is not obvious, because the macro model's tail does *not* track
  the real tail sweep by sweep (0 of 429 at r = 7) and only the top-cell parity
  agrees. But it is dominated by simply fitting the extender in lanes, so it is
  a curiosity rather than a deliverable.

## Status

Termination, correctness, and the physical -> Catalan reduction are solved on
paper but not formalized. *That* the parities differ is now settled — two
distinct finite-range slopes, `t = 9.0` (see Data) — but *why* they differ is
open, as is whether any subexponential algorithm exists. Why L = 69 is tiny is
diagnosed dynamically (13-cell phase resonance, phase 23 uniquely closing) but
not derived algebraically.

Odds of a subexponential algorithm existing and being findable from here:
**~1-2%**, essentially unmoved through every round. What has changed is that the
*structural* position is now much better than the algorithmic one. The project
now has a complete recursive lift `A_n -> A_{n+1}` with a seven-state local
carrier and `w < 3p`, an even-sector version with `w < 4p`, a two-line
reversible recursion for one whole sweep on plane trees, an explicit inverse, a
one-line theorem saying no merging quotient can ever work, and a height-
preserving reversor that turns the two standing reflection conjectures into
theorems. None of that has produced a program that beats simulating the orbit,
and three of the things that looked like they might have — a polynomial-size
itinerary grammar, sub-weight grammar evaluation, and `Xi`-weight multiplicity
counting — do not survive audit (see rules 6 and 7, and the dead-ends list).

The one route that keeps being pointed at is: an *injective* representation with
cheap random access, in the spirit of modular exponentiation rather than of a
smaller state. That is the only shape the no-go leaves open. Nothing measured so
far suggests one exists here: at L = 43 the projection of the canonical orbit to
`C_42` is injective on 2,232,313 of 2,232,314 phases, a `K`-step changes ~0.58N
recursive coordinates, and both the `Xi` and `Omega` recursions are 12-67x
*slower* than direct simulation when actually run.

Two of the open sub-problems this section used to list are now closed, both by
the pointed reversor `E_n` (Settled theory): the canonical `Xi`-weight
palindrome and the global height reflection
`delta_t = -delta_{p_L - 2L + 1 - t}` are theorems, not measurements, and
`eps_L` has an exact midpoint characterization. Neither buys any speed — `E_n`
gives the state at the opposite canonical phase without knowing `p_L`, but not
the *distance* between those phases, and at `O(L^2)` against `O(L)` for one
sweep it cannot even halve the walk. The unique-marker obstruction is untouched.

Still open, if anyone wants them: prove the 84-excursion palindrome at L = 69,
which is what the `43 = 62 + 110` decomposition rests on — the reversor proves
the canonical height and exceptional-ECO palindromes but does not yet act on the
`56|13` clean-cut defect classification; find a cheap `P -> P^dagger` prefix map
(without it, L = 69 is diagnosed but not explained); characterize which cuts
`n|s` carry the constant `n - s`, since it is neither all clean cuts nor exactly
the stem set; measure whether the residual tree sections reached from the
canonical chain form a `exp(o(n))` sub-language, which is the last unmeasured
thing the tree recursion asks for; and settle whether the two parity growth
*bases* really differ asymptotically, now that the finite-range slopes are known
to (see Data).

**Recommendation on record: treat the research thread as closed and take the
AVX2 hot loop as the deliverable.**
