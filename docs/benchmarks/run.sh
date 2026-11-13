#!/usr/bin/env bash
# Benchmark harness for docs/benchmark-plan.md. Runs inside WSL/Linux.
# Usage: RUNS=3 INPUTS="BER11_L_NOuser Motala3" SYNTH="1M 10M 50M" bash run.sh
# Resumable: a run whose .time log already exists is skipped.
set -euo pipefail

RUNS=${RUNS:-3}
BENCH=${BENCH:-$HOME/bench}
INPUTS=${INPUTS:-"BER11_L_NOuser Motala3"}
SYNTH=${SYNTH:-"1M 10M 50M"}
THREADS=${THREADS:-"1 2 4 8 12"}
OUT=${OUT:-$BENCH/results}
BIN=$BENCH/pmdshodh/target/release/pmdshodh
PMDTOOLS=$BENCH/pmdtools/PMDtools-0.60/pmdtools.0.60.py
DATA=$BENCH/data
WORK=$BENCH/work
export PATH=/root/miniforge3/envs/pmdshodh-bench/bin:$PATH

mkdir -p "$OUT/logs" "$WORK"
TIMES=$OUT/times.tsv
[ -s "$TIMES" ] || printf 'name\tinput\tconfig\trun\twall_s\tuser_s\tsys_s\tmax_rss_kb\texit\n' > "$TIMES"

log() { echo "[$(date -u +%FT%TZ)] $*" | tee -a "$OUT/progress.log" >&2; }

# Record the environment once.
if [ ! -s "$OUT/environment.txt" ]; then
  {
    echo "date_utc: $(date -u +%FT%TZ)"
    echo "kernel: $(uname -srm)"
    grep -m1 'model name' /proc/cpuinfo
    echo "nproc: $(nproc)"
    free -m | head -2
    echo "pmdshodh: $("$BIN" --version) commit $(git -C "$BENCH/pmdshodh" rev-parse HEAD)"
    echo "samtools: $(samtools --version | head -1)"
    echo "metaDMG: $(metaDMG-cpp pmd 2>&1 | grep -m1 version | sed 's/^[[:space:]]*-> //')"
    echo "metaDMG conda package: $(/root/miniforge3/bin/conda list -p /root/miniforge3/envs/pmdshodh-bench metadmg 2>/dev/null | awk '/^metadmg/{print $2}')"
    echo "python3: $(python3 --version)"
    echo "PMDtools archive md5: $(md5sum "$BENCH/pmdtools/0.60.tar.gz" | cut -d' ' -f1)"
    echo "runs per configuration: $RUNS"
    for i in $INPUTS; do echo "input $i: $(stat -c %s "$DATA/$i.bam") bytes md5 $(md5sum "$DATA/$i.bam" | cut -d' ' -f1)"; done
    echo "windows power scheme: $(powercfg.exe /getactivescheme 2>/dev/null | tr -d '\r' || echo unknown)"
  } > "$OUT/environment.txt"
fi

# timed NAME INPUT CONFIG RUN CMD...: run CMD under /usr/bin/time -v, append a row to times.tsv.
timed() {
  local name=$1 input=$2 config=$3 run=$4; shift 4
  local tlog="$OUT/logs/$name.time"
  if [ -s "$tlog" ]; then log "skip $name (done)"; return; fi
  log "start $name"
  printf '%s\n' "$*" > "$OUT/logs/$name.cmd"
  local rc=0
  /usr/bin/time -v -o "$tlog.tmp" bash -c "set -o pipefail; $*" 2> "$OUT/logs/$name.stderr" || rc=$?
  mv "$tlog.tmp" "$tlog"
  awk -v n="$name" -v i="$input" -v c="$config" -v r="$run" -v rc="$rc" -F': ' '
    /Elapsed \(wall clock\)/ { k=split($2,t,":"); w=(k==3)? t[1]*3600+t[2]*60+t[3] : t[1]*60+t[2] }
    /User time/ {u=$2} /System time/ {s=$2} /Maximum resident/ {m=$2}
    END { printf "%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\n", n,i,c,r,w,u,s,m,rc }' "$tlog" >> "$TIMES"
  log "done  $name rc=$rc"
  [ "$rc" -eq 0 ] || { log "FAILED $name; stopping"; exit 1; }
}

warm() { log "warm cache $1"; cat "$1" > /dev/null; }

bench_input() {
  local i=$1 in=$DATA/$1.bam
  warm "$in"
  for r in $(seq 1 "$RUNS"); do
    timed "$i.pmdtools_filter.r$r" "$i" pmdtools_filter "$r" \
      "samtools view -h $in | python3 -W ignore $PMDTOOLS --threshold 3 --header | samtools view -b -o $WORK/$i.pmdtools.bam -"
  done
  timed "$i.pmdtools_printDS.r1" "$i" pmdtools_printDS 1 \
    "samtools view -h $in | python3 -W ignore $PMDTOOLS --printDS > $WORK/$i.pmdtools.ds"
  for t in $THREADS; do for r in $(seq 1 "$RUNS"); do
    timed "$i.filter.t$t.r$r" "$i" "filter_t$t" "$r" \
      "$BIN filter --threshold 3 --quiet -@ $t $in -o $WORK/$i.filter.t$t.bam"
  done; done
  for t in 1 8; do for r in $(seq 1 "$RUNS"); do
    timed "$i.compat.t$t.r$r" "$i" "compat_t$t" "$r" \
      "$BIN filter --threshold 3 --quiet --pmdtools-compat -@ $t $in -o $WORK/$i.compat.t$t.bam"
    timed "$i.score.t$t.r$r" "$i" "score_t$t" "$r" \
      "$BIN score --quiet -@ $t $in -o $WORK/$i.score.t$t.tsv"
    timed "$i.metadmg.t$t.r$r" "$i" "metadmg_t$t" "$r" \
      "metaDMG-cpp pmd -@ $t $in > $WORK/$i.metadmg.t$t.txt"
  done; done
  timed "$i.score_compat.t8.r1" "$i" score_compat_t8 1 \
    "$BIN score --quiet --pmdtools-compat -@ 8 $in -o $WORK/$i.score_compat.t8.tsv"
}

# Synthetic inputs for the memory test: BER11_L_NOuser reads, repeated as needed.
make_synth() {
  local label=$1 n src=$DATA/BER11_L_NOuser.bam out=$WORK/synth_$1.bam
  case $label in 1M) n=1000000;; 10M) n=10000000;; 50M) n=50000000;; *) echo "bad size $label"; exit 1;; esac
  [ -s "$out" ] && return
  log "make $out ($n records)"
  { samtools view -H "$src"; for k in 1 2 3 4; do samtools view "$src"; done; } \
    | awk -v n="$n" '/^@/ {print; next} {print; if (++c >= n) exit}' \
    | samtools view -@ 4 -b -o "$out.tmp" - || true
  mv "$out.tmp" "$out"
  [ "$(samtools view -c "$out")" -eq "$n" ] || { log "synth $label has wrong record count"; exit 1; }
}

for i in $INPUTS; do bench_input "$i"; done
for s in $SYNTH; do
  make_synth "$s"; warm "$WORK/synth_$s.bam"
  for t in 1 8; do for r in $(seq 1 "$RUNS"); do
    timed "synth_$s.filter.t$t.r$r" "synth_$s" "filter_t$t" "$r" \
      "$BIN filter --threshold 3 --quiet -@ $t $WORK/synth_$s.bam -o $WORK/synth_$s.filter.t$t.bam"
  done; done
done
log "ALL DONE"
