#!/usr/bin/env bash
set -euo pipefail

# First-run driver for the compact X1(37), F_2 degree-9 translate test.
#
# Usage:
#   ./run_first.sh [CORES] [NCHUNKS]
#
# The driver verifies the hard-coded map and the cached order-5 generator once,
# then runs every degree-type chunk.  Worker shortcuts are used only after the
# verification log contains the expected success markers.

CORES="${1:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)}"
NCHUNKS="${2:-16}"
MAIN="${MAIN:-X1_37_prop49_F2_short_complete.m}"
LOGDIR="${LOGDIR:-logs_first_run}"
VERBOSE_EVERY="${VERBOSE_EVERY:-1000}"
MAX_SURVIVORS="${MAX_SURVIVORS:-1}"
GEN_I="${GEN_I:-1}"
GEN_J="${GEN_J:-4}"

if ! command -v magma >/dev/null 2>&1; then
    echo "ERROR: magma was not found on PATH." >&2
    exit 2
fi

if [[ ! -s X1_37_X0plus_map_data.m ]]; then
    echo "ERROR: X1_37_X0plus_map_data.m is missing or empty." >&2
    exit 2
fi

mkdir -p "$LOGDIR"

echo "======================================================================"
echo "One-time verification"
echo "======================================================================"
magma -b RUN_TRANSLATE_TEST:=false \
         VERIFY_MAP_BY_RECONSTRUCTION:=true \
         VERIFY_MODEL_BASIC:=true \
         VERIFY_DELTA_ORDER:=true \
         F2_CACHED_SCAN_TARGETS_ONLY:=false \
         CACHED_GENERATOR_I:="$GEN_I" \
         CACHED_GENERATOR_J:="$GEN_J" \
         "$MAIN" 2>&1 | tee "$LOGDIR/verify.log"

grep -q "MAP_VERIFICATION: SUCCESS" "$LOGDIR/verify.log" || {
    echo "ERROR: map verification did not report success." >&2
    exit 3
}

grep -q "cyclic_delta_ok[[:space:]]*=[[:space:]]*true" "$LOGDIR/verify.log" || {
    echo "ERROR: order-5 quotient delta was not verified." >&2
    exit 3
}

grep -q "Using candidate fiber pair <${GEN_I},${GEN_J}> as a verified order-5 generator" "$LOGDIR/verify.log" || {
    echo "ERROR: cached worker generator <${GEN_I},${GEN_J}> was not verified." >&2
    exit 3
}

echo "======================================================================"
echo "Running degree-type chunks"
echo "CORES=$CORES"
echo "NCHUNKS=$NCHUNKS"
echo "LOGDIR=$LOGDIR"
echo "======================================================================"

seq 0 $((NCHUNKS - 1)) | xargs -P "$CORES" -I {} bash -c '
    idx="$1"
    logdir="$2"
    nchunks="$3"
    main="$4"
    verbose_every="$5"
    max_survivors="$6"
    gen_i="$7"
    gen_j="$8"
    logfile="$logdir/chunk_${idx}.log"

    echo "Starting chunk $idx at $(date)" > "$logfile"
    set +e
    magma -b RUN_TRANSLATE_TEST:=true \
             VERIFY_MAP_BY_RECONSTRUCTION:=false \
             VERIFY_MODEL_BASIC:=false \
             VERIFY_DELTA_ORDER:=false \
             QUIET_SETUP:=true \
             F2_CACHED_SCAN_TARGETS_ONLY:=true \
             CACHED_GENERATOR_I:="$gen_i" \
             CACHED_GENERATOR_J:="$gen_j" \
             NCHUNKS:="$nchunks" \
             CHUNK_INDEX:="$idx" \
             VerboseEvery:="$verbose_every" \
             MaxSurvivors:="$max_survivors" \
             "$main" >> "$logfile" 2>&1
    status=$?
    echo "Finished chunk $idx at $(date) with status $status" >> "$logfile"
    exit "$status"
' _ {} "$LOGDIR" "$NCHUNKS" "$MAIN" "$VERBOSE_EVERY" "$MAX_SURVIVORS" "$GEN_I" "$GEN_J"

bad=0
for idx in $(seq 0 $((NCHUNKS - 1))); do
    logfile="$LOGDIR/chunk_${idx}.log"
    if grep -q "CHUNK_RESULT: SUCCESS" "$logfile"; then
        echo "chunk $idx: SUCCESS"
    else
        echo "chunk $idx: NOT SUCCESS -- see $logfile"
        bad=1
    fi
done

if [[ "$bad" -eq 0 ]]; then
    echo "GLOBAL_RESULT: SUCCESS"
else
    echo "GLOBAL_RESULT: INCONCLUSIVE_OR_FAILED"
    exit 4
fi
