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
`{1:56, 3:28, 5:14, 9:8, 11:4, 13:1}` summing to 339. Only position 111 (the
"hole") varies, so `p_L = 451 * (q_{L-1} / 112) + sum(holes)`. This gives
weights and never states, which is why it does not speed anything up.

## Data

`B_L` for L = 4..46:

    3, 9, 9, 22, 23, 23, 23, 54, 23, 117, 185, 211, 451, 451, 451, 918, 451,
    1854, 919, 1853, 919, 3729, 4667, 5608, 7493, 11251, 4673, 115288, 6559,
    36504, 8887, 201627, 37551, 386120, 15927, 705934, 21567, 382512, 185403,
    2232272, 38513, 729845, 33767

Plus `B_60 = 7,640,649`, `B_64 = 22,789,705`, `B_66 = 74,827,985`,
**`B_69 = 1,124,943`**, and `p_49 = 94,695,837`.

**`B_L` is not smooth.** Almost every value has a prime factor above
`sqrt(B_L)`; `B_15 = 211`, `B_22 = B_24 = 919`, `B_29 = 11251`, `B_30 = 4673`,
`B_34 = 8887` are prime, and `B_69 = 3 * 374981` with 374981 prime. So `B_L` is
*not* an lcm of small independent periods, which kills the shape of every
decomposition route.

**Parity separation is robust; two distinct growth bases are not.** Odd `B_L`
dwarfs even at the same `L` (`B_43 = 2,232,272` vs `B_44 = 38,513`, 58x). But
least-squares on `log B_L` gives odd 1.3575 (L >= 16) / 1.3553 (L >= 24) against
even 1.2004 / 1.1762 — while adding `B_60/64/66` pushes the even base to
**1.3471** over L = 40..66. The fit moves by ~0.15 with the window, and a common
base with different constants is not excluded by data to L = 46. **Distrust any
quoted pair** (1.37/1.23, 1.42/1.23, ...); short-window fits at L <= 44 keep
re-endorsing them and that is a documented artifact. Mean checkpoint gap is ~4.1
for both parities, so odd extenders are not slow per checkpoint — they visit
more checkpoints.

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
  8, 112, swept to 45. The law `H^112(1^16 u) = 1^16 J(u)` is real, but `J` has
  no closed form.
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
- **Global palindrome involution `J A J = A^{-1}`** — vacuous as stated, since
  every permutation is strongly real; and preserving every prefix-reset set needs
  one common centre, which fails (at L = 16, `S_L` forces c = 0 but `S_1` gives
  450, and `S_9` has no centre at all). The palindromes are real but live in
  *induced* time, one involution per level: per-level centres at L = 16,
  k = 1..8 are 450, 448, 448, 442, 428, 428, 428, 428, with gaps `2^j - 2`.
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
   renormalized map `R_r`. Measured flat-vs-checkpoint gate counts at
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

## Status

Termination, correctness, and the physical -> Catalan reduction are solved on
paper but not formalized. Why the parities differ, and whether any
subexponential algorithm exists, are open. Why L = 69 is tiny is diagnosed
dynamically (13-cell phase resonance, phase 23 uniquely closing) but not derived
algebraically.

Odds of a subexponential algorithm existing and being findable from here: **~2%**,
with one route left — derive the growth bases, which needs an idea rather than a
measurement. Open sub-problems, if anyone wants them: prove the phase-reflection
law abstractly and find a cheap `P -> P^dagger` prefix map (without it, L = 69 is
diagnosed but not explained); and for the general fast-forward, the only
surviving target named is a nonlinear defect-cocycle representation closed under
composition.

**Recommendation on record: treat the research thread as closed and take the
AVX2 hot loop as the deliverable.**
