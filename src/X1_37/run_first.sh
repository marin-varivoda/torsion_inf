#!/usr/bin/env bash
set -euo pipefail

CORES="${1:-16}"
NCHUNKS="${2:-16}"
LOGDIR="${LOGDIR:-logs_minimal_run}"
mkdir -p "$LOGDIR"

command -v magma >/dev/null 2>&1 || { echo "ERROR: magma not found" >&2; exit 2; }
[[ -s X1_37_X0plus_map_data.m ]] || { echo "ERROR: missing X1_37_X0plus_map_data.m" >&2; exit 2; }

echo "======================================================================"
echo "One-time verification"
echo "======================================================================"
magma verify.m 2>&1 | tee "$LOGDIR/verify.log"

grep -q "MAP_VERIFICATION: SUCCESS" "$LOGDIR/verify.log"
grep -q "ST_RELATION_VERIFICATION: SUCCESS" "$LOGDIR/verify.log"
grep -q "DELTA_ORDER: SUCCESS" "$LOGDIR/verify.log"

echo "======================================================================"
echo "Running chunks: CORES=$CORES NCHUNKS=$NCHUNKS LOGDIR=$LOGDIR"
echo "======================================================================"

seq 0 $((NCHUNKS - 1)) | xargs -P "$CORES" -I {} bash -c '
    idx="$1"; logdir="$2"; nchunks="$3"
    logfile="$logdir/chunk_${idx}.log"
    echo "Starting chunk $idx at $(date)" > "$logfile"
    set +e
    magma -b NCHUNKS:="$nchunks" CHUNK_INDEX:="$idx" worker.m >> "$logfile" 2>&1
    status=$?
    echo "Finished chunk $idx at $(date) with status $status" >> "$logfile"
    exit "$status"
' _ {} "$LOGDIR" "$NCHUNKS"

bad=0
for idx in $(seq 0 $((NCHUNKS - 1))); do
    if grep -q "CHUNK_RESULT: SUCCESS" "$LOGDIR/chunk_${idx}.log"; then
        echo "chunk $idx: SUCCESS"
    else
        echo "chunk $idx: NOT SUCCESS -- see $LOGDIR/chunk_${idx}.log"
        bad=1
    fi
done

if [[ "$bad" -eq 0 ]]; then
    echo "GLOBAL_RESULT: SUCCESS"
else
    echo "GLOBAL_RESULT: INCONCLUSIVE_OR_FAILED"
    exit 4
fi
