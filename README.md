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

The theory behind it, the exact values of `B_L` up to L = 69, and — more useful
to anyone picking this up — the list of approaches that were tried and measured
dead, are written up in [`docs/research-notes.md`](docs/research-notes.md). The
programs that check each claim are in [`research/`](research/README.md). The
short version is that no subexponential algorithm was found and the estimated
odds of one being reachable from here are about 2%.

## Credit

Development into Snaperz extenders spans multiple years. If you played a role but we forgot to list you here, contact us and we will add your name.

- Snaperz
- Space Walker
- MadCloud101
- G4me4u
- Ralp
