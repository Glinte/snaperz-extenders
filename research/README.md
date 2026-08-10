# Research programs

One-off programs written while working out exact pulse counts for period-12
Snaperz extenders. They are not part of the simulator and nothing in `src/`
depends on them; they exist so that every claim in
[`../docs/research-notes.md`](../docs/research-notes.md) can be re-derived from
scratch rather than taken on trust.

Read the notes first. This file only says which program checks what.

## Building and running

The C programs are built by the `Makefile` here, not by the top-level CMake
build, and `cmake` is not required:

```bash
make            # everything, into ./build
make build/r14a # just one
make clean
```

The Python scripts need only a Python 3.12 stdlib. `snaperz_p12.py` will
optionally shell out to a C compiler for its hot loop and falls back to pure
Python if none is found. Scripts under `audit/` import `audit/core.py`, so run
them from their own directory:

```bash
python3 snaperz_p12.py --self-test
cd audit && python3 r14_small.py
```

Generated outputs (`build/`, `sections69.txt`, `g18.json`) are ignored by git;
everything regenerates in seconds to a couple of minutes.

## `snaperz_p12.py` — the reference implementation

The whole derivation in one file: Catalan states, the blockless sweep `A_L`, the
good-pair fibre argument giving `T_L = (1 + eps_L) p_L - (L - 1)`, the
checkpoint map `H`, and the constant-return-time induced sections.

```bash
python3 snaperz_p12.py 30 31      # T_L for a range of lengths
python3 snaperz_p12.py --table 40 # the table of B_L, p_L, eps_L
python3 snaperz_p12.py --self-test
python3 snaperz_p12.py --no-c     # force the pure-Python path
```

`--self-test` is the one that matters: it checks `T_L` against a direct
simulation of the physical extender for L = 1..22, checks the C core against the
Python fallback for L = 1..24, and confirms the constant-time sections
`{1: 1, 3: 2, 7: 8, 8: 8, 16: 112}`.

`checkpoint_step` is the `H` map and `orbit_data` returns `(q, p, eps)`.

## `cumulative/` — the 3-state carry automaton

The rewrite onto cumulative coordinates `s_i = sum_{j <= i} a_j`, which is what
made `src/snaperz_extender_fsm.h` possible.

- `snaperz_cumsweep.c` — four cores over the same orbit (a-coordinate branchy,
  a-coordinate branchless, cumulative branchless, cumulative + FSM + SSE),
  cross-checked against each other and timed. `./build/snaperz_cumsweep 49`
  reproduces the 35.1s / 27.9s / 7.25s figures in the notes.
- `snaperz_cumsweep_verify.py` — checks the cumulative sweep against
  `snaperz_p12.py` on random Catalan states, and reports the `a_i == 0` fraction
  along the orbit, which is why the branchless variant alone does not pay.

Note the measured conclusion: the cumulative rewrite on its own is a
*pessimisation*. The win is the factorization it enables.

## `diagonals/` — is `A^t` compressible?

To evaluate `A^t` in one left-to-right pass you must carry the anti-diagonal, so
the number of distinct anti-diagonals `N(t)` lower-bounds any transducer for
`A^t`. Polynomial `N(t)` would give `B_L` by repeated squaring; exponential
closes every fast-forward route that composes `A` with itself.

- `snaperz_diag.c` — sampled count along one orbit prefix, HyperLogLog so memory
  is constant. `./build/snaperz_diag <length> <sweeps> <max_t>`
- `snaperz_diag_exact.c` — no sampling: enumerates the whole state space, so
  every window of `t` consecutive rows is counted. This is the one to trust; the
  sampled version bent the growth curve downward three separate times (see rule
  3 in the notes). `./build/snaperz_diag_exact <length>`

Answer: `N(t)` grows at least `t^4` with the exponent rising with `L`, and
`N(t) <= B_L * L` means the transducer is the size of the orbit it would
replace. Route closed.

## `audit/` — independent reconstruction and the round-13/14 audits

`core.py` rebuilds `A_L` from the gate rule alone, with no shared code with
`snaperz_p12.py`. Everything else here is built on that reconstruction, so the
two implementations agreeing is a real cross-check. The C programs each inline
the same six-line `Astep` for the same reason.

Round 13 — auditing an external write-up:

- `big.c` — `B_L` over a range, used to validate the reconstruction against the
  data table in the notes (L = 4..44 match exactly).
- `s3.py` — the decorated parity quotient as an invariant subsystem: good pairs,
  the ternary-tree bijection, and the ordered-pair counts `N_n`.
- `s56.py` — the checkpoint map `H` on plane trees, `z(H(b)) = b_0`, and the
  `Theta`/`G` factorization.
- `probe.c` — is `ell >= 56` a property of the 16-clock section only, or of every
  checkpoint? (Only the section; off-section `min ell = 0`.)
- `where.c` — locates the 16 sub-56 checkpoints, which turn out to be the escape
  that *is* the closure.
- `extra.c` — traces `ell` across the claimed reflection at outer step 780/781,
  for L = 69 and its odd neighbours.
- `rank.c` — the two affine ranks over `F_2`. These reproduce, but the inference
  drawn between them in the external write-up does not hold; see the notes.

Round 14 — auditing the continuation:

- `r14_small.py` — the small combinatorial claims: `z(H(b)) = b_0`, `G`
  preserving root degree, the subtree-size multiset failing, `p_13 = 129`, the
  15 `Theta`-compatible phases, and `p_18 = 468` with the 112-gap word.
- `r14a.c` — the checkpoint-time reflection identity
  `t_k + t_{N-k} = p_L - L + z(x_k)` (0 failures over 268,129 indices at
  L = 69), the 56+13 clean-section census, the phase-reflection law, and the
  universal terminal cone. Writes `sections69.txt`.
- `r14b.c` — the lazy-winching reformulation, defect density, rotating-frame
  states, branch vectors, and the linear complexity of `(A^t E)_0 mod 2`. Four
  of the five algorithmic eliminations are here.
- `r14c.c` — pins down the `L_j = W_j . S_j` factorization and redoes the
  branch-vector count with the branch read *during* the sweep.
- `r14d.c` — the repaired factorization (swap only when the ceiling is slack)
  and the L = 69 separator word's palindrome decomposition.
- `r14e.c` — the controlled 129-phase experiment: only phase 23 reaches `E_69`,
  at exactly return 1,614. The strongest L = 69 result.
- `r14f.c` — the `56 + s` tail phase-lock table and the injected root `q`, which
  is where the 16-barrier shows up again at `s = 17`.

Round 15 — auditing the `R_r` renormalization:

- `r15a.py` — the state-level renormalization itself: the closed form for `H`,
  `R_r = H_r . S_{30,r} . A_r^435` in both spellings, the literal 112-return and
  `466 + c` full-sweep identities, the first failure of the naive
  `H_r . A_r^435` at r = 6, the periods `M_r` and the recurrence
  `p_{16+r} = 466 M_r + sum z(y_j)` against direct simulation, the
  `Psi = H.Theta` reflection candidate failing at m = 6, and the section
  `Gamma_s`.
- `r15b.c` — the reflection constant over every clean `n|s` cut of every orbit it
  can run. This is what shows the proposed general law is false at generic cuts
  and true at the stem cuts, and that `43 = n - s` with `B_s` cancelling.
  `./build/r15b <L> <n>` tests a single cut.
- `r15c.c` — the lock offset `(phi_k - t_k) mod p_s` over every clean checkpoint,
  which kills the phase-lock proof of the reflection law: 62 distinct offsets at
  L = 69, and the constant it would force is 93 rather than 43.
- `r15d.py` — welds `R_r` to the 112-block template: block count `= M_r`, hole
  `= 16 + z(y_j)`, universal 111-gap prefix summing to 450.

- `r15e.py` — the recursive-16-peel test, which is the one route the `R_r`
  derivation suggested and the one it kills: the `R_r` orbit never returns to
  shielded form, so there is no second level to peel.
- `r15f.py` — the cycle-quotient audit for `Gamma_s`, which is the decisive test
  for phase coordinates. `Gamma_s` scatters the canonical `A_s`-cycle across
  ~`p_s` distinct target cycles, so the cycle identifier is as large as the
  state.
- `r15g_bench.cpp` — does `R_r` pay? It checks the renormalized counter against
  `blockless::run_scalar` on `B_L`, `eps_L` and `T_L`, times the two, and
  measures the wavefront's ns/sweep across L. Answers: yes on the scalar path
  (1.2-3.1x), no on the wavefront, whose cost does not depend on L at all.
  `./build/r15g_bench renorm <L>` runs the renormalized counter alone.
- `r15h_avx512.cpp` — validates `src/blockless_avx512.h`, the 64-lane windows
  that carry L <= 127, against the 256-bit windows on every length where both
  work and against the `B_L` table.

The two C++ programs include `../../src/blockless.h`, so they check the shipped
counter rather than a copy of it.

## What is not here

The Lean formalization lives in a separate repository. Its early attempts
formalized a wrong, simplified transition and do not count toward anything
claimed in the notes.
