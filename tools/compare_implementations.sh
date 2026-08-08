#!/usr/bin/env bash
#
# Builds every extender implementation at a range of lengths and checks that
# they all report the same pulse count, then times them at one longer length.
#
#   ./tools/compare_implementations.sh [max_verify_length] [bench_length]
#
# Implementations that the host CPU cannot run are skipped.

set -u

BASEDIR=$(cd "$(dirname "$0")/.." && pwd)
MAX_VERIFY=${1:-40}
BENCH=${2:-49}
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

CXXFLAGS="-O3 -march=native -std=c++17 -DCHECK_LOOP=1 -DLOG_STATUS_UPDATES=0"

IMPLS="fallback fsm"
if gcc -march=native -dM -E - </dev/null | grep -q "__AVX2__"; then
  IMPLS="$IMPLS avx2"
else
  echo "note: no AVX2 on this CPU, skipping the AVX2 implementation"
fi

flags_for()
{
  case "$1" in
    fallback) echo "-DSNAPERZ_FORCE_FALLBACK=1" ;;
    fsm)      echo "-DSNAPERZ_FORCE_FSM=1" ;;
    avx2)     echo "" ;;
  esac
}

echo "Verifying lengths 2..$MAX_VERIFY"
for impl in $IMPLS; do
  : > "$OUT/$impl"
  for n in $(seq 2 "$MAX_VERIFY"); do
    # shellcheck disable=SC2086
    if ! g++ $CXXFLAGS -DSNAPERZ_LENGTH="$n" $(flags_for "$impl") \
         -o "$OUT/bin_$impl" "$BASEDIR/src/main.cpp" 2>/dev/null; then
      echo "$n build-failed" >> "$OUT/$impl"
      continue
    fi
    echo "$n $("$OUT/bin_$impl" | tail -1 | awk '{print $2}')" >> "$OUT/$impl"
  done
done

status=0
reference=$(echo "$IMPLS" | awk '{print $1}')
for impl in $IMPLS; do
  [ "$impl" = "$reference" ] && continue
  if diff -q "$OUT/$reference" "$OUT/$impl" >/dev/null; then
    echo "  $impl agrees with $reference"
  else
    echo "  $impl DISAGREES with $reference:"
    diff "$OUT/$reference" "$OUT/$impl" | head -20
    status=1
  fi
done

echo
echo "Timing length $BENCH"
for impl in $IMPLS; do
  # shellcheck disable=SC2086
  g++ -O3 -march=native -std=c++17 -DCHECK_LOOP=0 -DLOG_STATUS_UPDATES=0 \
      -DSNAPERZ_LENGTH="$BENCH" $(flags_for "$impl") \
      -o "$OUT/bin_$impl" "$BASEDIR/src/main.cpp" 2>/dev/null || continue
  start=$(date +%s.%N)
  pulses=$("$OUT/bin_$impl" | tail -1 | awk '{print $2}')
  end=$(date +%s.%N)
  elapsed=$(awk -v a="$start" -v b="$end" 'BEGIN { printf "%.2f", b - a }')
  printf "  %-9s %12s pulses  %8s s\n" "$impl" "$pulses" "$elapsed"
done

exit $status
