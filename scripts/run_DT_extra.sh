#!/bin/sh

root=$(cd "$(dirname "$0")/.." && pwd)
src=$root/src
out=$root/out/DT/extra_X1_2_N

mkdir -p "$out"
cd "$src" || exit 1

if ! command -v magma >/dev/null 2>&1; then
	echo "magma not found in PATH" >&2
	exit 1
fi

pids=
total=0

for name in \
	gonality_X1_2_14 \
	gonality_X1_2_16 \
	gonality_X1_2_18 \
	gonality_X1_2_20
do
	echo "starting $name"
	magma -n "DT/extra_X1_2_N/${name}.m" >"$out/${name}.out" 2>&1 &
	pids="$pids $!"
	total=$((total + 1))
done

echo "started $total jobs"

tick=0
while :; do
	running=0
	for pid in $pids; do
		if kill -0 "$pid" 2>/dev/null; then
			running=$((running + 1))
		fi
	done
	done=$((total - running))
	case $((tick % 4)) in
	0) spin='|' ;;
	1) spin='/' ;;
	2) spin='-' ;;
	3) spin='\' ;;
	esac
	tick=$((tick + 1))
	printf '\rProgress: %d of %d  %s' "$done" "$total" "$spin"
	if [ "$running" -eq 0 ]; then
		break
	fi
	sleep 1
done
printf '\rProgress: %d of %d\n' "$total" "$total"

failed=0
for pid in $pids; do
	if ! wait "$pid"; then
		failed=$((failed + 1))
	fi
done

if [ "$failed" -gt 0 ]; then
	echo "$failed job(s) failed" >&2
	exit 1
fi

echo "all done"
