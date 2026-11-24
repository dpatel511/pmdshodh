#!/usr/bin/env bash
# One condition of the damage-curve study (docs/damage-curve-study.md): simulate, score with
# every scorer, evaluate. Writes NAME.* in the current directory; NAME.eval.json is the result.
#
#   PMDSHODH=/path/to/pmdshodh tests/sim/run_case.sh base --seed 11 --divergence 0.001
#
# Arguments after NAME go to simulate.py. pmdshodh gets the same --library and --treatment.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
pm=${PMDSHODH:-pmdshodh}
py=${PYTHON:-python3}
name=$1
shift

$py -I "$here/simulate.py" --out "$name" "$@"
read -r lib treat < <($py -I -c 'import json,sys; p=json.load(open(sys.argv[1])); print(p["library"], p["treatment"])' "$name.params.json")
mode=(--library "$lib" --treatment "$treat")

"$pm" damage --quiet "${mode[@]}" "$name.sam" -o "$name.damage.json"
$py -I "$here/fit.py" "$name.damage.json" > "$name.fit.json"
q=$($py -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["q"])' "$name.fit.json")

"$pm" score --quiet "${mode[@]}" "$name.sam" -o "$name.default.tsv"
"$pm" score --quiet "${mode[@]}" --PMDpparam "$q" "$name.sam" -o "$name.fitted.tsv"
$py -I "$here/score.py" "$name" --curve true > "$name.true.tsv"
$py -I "$here/score.py" "$name" --curve true --null exact > "$name.oracle.tsv"
$py -I "$here/score.py" "$name" --curve true --null exact --length > "$name.oracle_length.tsv"

$py -I "$here/evaluate.py" \
    default="$name.default.tsv" fitted="$name.fitted.tsv" true="$name.true.tsv" \
    oracle="$name.oracle.tsv" oracle_length="$name.oracle_length.tsv" > "$name.eval.json"
