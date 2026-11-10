# Simulating the Snaperz Piston Extender
A simple C++ program for simulating a Snaperz Piston Extender, a specific subclass of piston extenders in Minecraft.

## Building the project
In order to build the project, run the `rebuild.sh` script.
This script contains the following commands, which will build the project in release mode:
```bash
mkdir -p "./build"
cd "./build"
cmake -DCMAKE_BUILD_TYPE=Release ..
make all
```
The project requires that `GCC` supporting C++17, `cmake`, and `make` are installed on your system.

In order to compile in debug mode, replace the third command with:
```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

### Linux
You likely already have the required build tools, simply run the commands above.

### Windows
If you are using Windows, the build dependencies can be installed through e.g. `cygwin64`, where the above command is run through the cygwin terminal.

## Running the project
After building the project, the simulation can be run through the following terminal command:
```bash
./build/extender
```
Depending on the size of the extender, this could take a significant amount of time. Be patient!

## Blazingly fast AVX2
The simulation also has AVX2 support, developed by G4me4u. This makes the program slightly less simple but at the same time blazingly fast! This feature requires AVX2 support on your CPU, and will otherwise use the traditional fallback implementation. CPU support is checked by running the command below in the terminal.
```bash
gcc -mavx2 -dM -E - < /dev/null | egrep "SSE|AVX" | sort
```
Your CPU supports AVX2 if the `#define __AVX2__ 1` line is shown in the output.

## SSE2 for older machines
Without AVX2 the program used to drop all the way down to the scalar fallback. There is now an intermediate implementation that needs only SSE2, which every x86-64 chip has. It stores prefix sums of the segment lengths rather than the lengths themselves, which turns a pulse into a three state finite automaton over the segments and lets most of the work be vectorized. The comment at the top of `src/snaperz_extender_fsm.h` has the derivation.

It is selected automatically when AVX2 is unavailable, and only for period 12, where the push limit is one. On a length 49 extender it takes about 15s against the fallback's 22s. It is not a replacement for the AVX2 implementation, which does the same run in 1.4s: that one parallelizes across concurrent pulses, which is a fundamentally better decomposition than parallelizing within a single pulse.

To check the implementations against each other and time them on your own machine:
```bash
./tools/compare_implementations.sh 44 49
```
A specific implementation can be forced with `-DSNAPERZ_FORCE_FSM=1` or `-DSNAPERZ_FORCE_FALLBACK=1`, and the extender configured from the build with `-DSNAPERZ_LENGTH=`, `-DSNAPERZ_PERIOD=` and `-DSNAPERZ_HARD_PUSH_LIMIT=`.

## Counting pulses without simulating them

For period 12 the pulse count can be derived rather than simulated, at least in
part. `src/blockless.h` and the `blockless` target implement the good-pair
reduction: deleting the extender block leaves a sweep on Catalan states, and the
pulse count comes back as `T_L = (1 + eps_L) * p_L - (L - 1)` from the cycle
length of that sweep. It is a large constant-factor win over simulating pulses,
not an asymptotic one, and it agrees with the AVX2 implementation everywhere
both have been run.

Its wavefront keeps the whole extender in vector lanes, so what matters is
whether the extender fits, not how long it is: the cost is about 10ns per sweep
whether L is 21 or 64. With 256-bit windows it fits L up to 64, and anything
longer dropped to the scalar sweep, roughly 20x slower. `src/blockless_avx512.h`
widens the windows to 64 lanes where AVX512BW is available, carrying L up to
127 — at L = 66 that is 1.3s against 25s. It is used only for lengths the narrow
windows cannot hold, since it runs at about 0.7x their speed where both apply.

The theory behind it, the exact values of `B_L` up to L = 69, and — more useful
to anyone picking this up — the list of approaches that were tried and measured
dead, are written up in [`docs/research-notes.md`](docs/research-notes.md). The
programs that check each claim are in [`research/`](research/README.md). The
short version is that no subexponential algorithm was found and the estimated
odds of one being reachable from here are about 2%.

## Half-cycle counting (period 12)

`./build/blockless --half` opts into the exact half-cycle counter for L > 16.
It uses the certified 16-coordinate reduction, stops at a reflection-fixed
edge or vertex, reconstructs the weighted period with its endpoint correction,
and obtains the physical parity by a short inverse rewind.

The reduced counter now vectorizes the 435-sweep blocks, tail sweeps, head
injections, particle-retiming reflection and midpoint rewind. Its finite
wavefront drains to a synchronous state at each checkpoint. For reflection,
suffix sums travel along the inverse wavefront and select each particle at its
required inverse time, avoiding a sequential reconstruction of all those rows.
The scalar implementation remains available as `half::run<false>(L)` for checks.

Measured end-to-end median seconds on an AMD Ryzen 5 PRO 215, compiled with
`g++ -O3 -march=native -std=c++17`. Each row uses five runs pinned to CPU 2, alternating order,
and checks the complete result against the existing counter:

| Length | SIMD half-cycle | Existing SIMD | Speedup |
| --- | ---: | ---: | ---: |
| 49 | 0.391 | 0.589 | 1.50x |
| 50 | 0.00173 | 0.00249 | 1.44x |
| 60 | 0.0331 | 0.0438 | 1.32x |
| 64 | 0.104 | 0.135 | 1.30x |
| 66 | 0.338 | 0.684 | 2.02x |
| 69 | 0.00507 | 0.0103 | 2.04x |

The existing counter uses AVX2 through 64 and AVX512 at 66 and 69. These are
machine-specific measurements on a shared, busy machine. An unpinned rerun
ranged from 0.97x at length 64 to 1.67x at 69; scheduling and frequency effects
matter, especially for short runs. The original scalar
half-cycle prototype took about 3.24s at length 66. SIMD makes the theoretical
saving useful against the already optimized baseline, although pipeline and
reflection overhead prevent a universal 2x speedup.

AVX2 handles up to 64 coordinates per reduced sweep; AVX512 extends that to
128 where available. Both require mass <= 127 (including injections) to keep
signed-byte arithmetic valid; larger inputs fall back to scalar operations.
The reduced orbit has width and mass L-16, so all its operations fit AVX2
through L=80, and AVX512 through L=143. The decoded physical midpoint can need
wider registers or scalar fallback sooner. Wide-kernel correctness is checked
at the boundaries; complete orbit timings above cover only the listed lengths.
The default counter remains unchanged; select this implementation with `--half`.

Reproduce the checks and comparisons with:

```bash
make -C research build/half_simd_test build/half_bench
./research/build/half_simd_test  # arbitrary states, boundaries and scalar certificates
./research/build/half_bench     # reference checks 17..50 and scalar/SIMD timing
taskset -c 2 ./research/build/half_bench simd # five-run medians on an allowed CPU
./research/build/half_bench 66   # all three counters at length 66
```

Tests cover forward/inverse blocks, injections, 16,641 reflection cases,
register and mass boundaries, both meeting types, odd/even periods and nonzero
epsilon. Native AVX512 and AVX2-only builds pass; address/undefined-behavior
sanitizers pass with leak detection disabled for the tracing environment.
The length-50 certificate remains 468 family steps, period 439822, rewind 15
and 439773 physical pulses.

The implementation follows the [half-cycle derivation](https://chatgpt.com/space/page_897a9c1bb6448191a5c9e9ebc112ff6a)
and `family_reversal.py` from the accompanying subexponential research.
This counts completion pulses; it does not produce an intermediate physical
configuration for every pulse, and still traverses half the reduced orbit.

## Credit

Development into Snaperz extenders spans multiple years. If you played a role but we forgot to list you here, contact us and we will add your name.

- Snaperz
- Space Walker
- MadCloud101
- G4me4u
- Ralp
