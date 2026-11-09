#!/bin/bash
# Runs an A/B experiment inside a Slurm job on one node; mkjob.sh submits it. Builds each variant's `measure` (or
# `DRIVER`) like CMake does, JCC padding included.
# Usage: ab-run.sh <job-dir>
#
# <job-dir> holds:
#   variants: lines of `name source [extra compile flags...]`; the source is ~/kh-perf/src/<source>/.
#   configs:  lines of `label policy threads driver-args...`, with policy one of local, interleave, all, preferred,
#             none; driver arguments come after `--threads` and `--reps`, so they can override the repetitions.
#   params:   shell assignments, e.g. ROUNDS=2 REPS=3 BUDGET_MIN=52 DRIVER=measure CPU_MODEL="Platinum 8276"; the
#             job exits unless `lscpu`'s model name contains CPU_MODEL.
# Outputs go to <job-dir>/out/<label>/<variant>.r<round>.tsv; runs with an output are skipped, so a resubmission
# continues where an earlier job stopped.

set -u
jd="$1"
ROUNDS=1; REPS=1; BUDGET_MIN=52; DRIVER=measure; CPU_MODEL="Platinum 8276"
source "${jd}/params"
t0=$(date +%s)

echo "== node"; hostname
lscpu | grep "^Model name" | grep -q "${CPU_MODEL}" || { echo "Not a ${CPU_MODEL} node: $(lscpu | grep '^Model name')"; exit 3; }
lscpu | grep -E "^(Model name|Socket|Core|Thread|NUMA node)"; echo "THP $(cat /sys/kernel/mm/transparent_hugepage/enabled) defrag $(cat /sys/kernel/mm/transparent_hugepage/defrag)"
numactl -H | grep -E "size|free"; grep -E "^(MemAvailable|Cached)" /proc/meminfo; echo "numa_balancing $(cat /proc/sys/kernel/numa_balancing)"
echo "$(g++ --version | head -1); -march=native: $(g++ -march=native -Q --help=target | grep -E '^  -march=' | awk '{ print $2 }')"
s0=$(lscpu -p=NODE,SOCKET | grep -v "^#" | awk -F, '$2 == 0 { print $1 }' | sort -un | paste -sd,)
echo "socket 0 NUMA nodes: ${s0}"

echo "== build"
b=/tmp/kh-${SLURM_JOB_ID:-$$}; mkdir -p ${b}
while read -r name src flags; do
    [ -z "${name}" ] && continue
    rsync -a --exclude "build*/" --exclude external/ ~/kh-perf/src/${src}/ ${b}/${name}/
    # The definitions are older trees': the headers set their own since v0.8.0.
    ( g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -funroll-loops -march=native -Wa,-mbranches-within-32B-boundaries \
        -DXXH_INLINE_ALL -DXXH_STATIC_LINKING_ONLY -DXXH_FORCE_MEMORY_ACCESS=0 -DDEVELOP -DL1_CACHE_LINE_SIZE=64 \
        ${flags} -I ${b}/${name}/include -I ${b}/${name}/support -I ${b}/${name}/benchmark -I ${b}/${name}/tools ${b}/${name}/benchmark/${DRIVER}.cpp -o ${b}/${name}.bin \
        ~/kh-perf/jemalloc/libjemalloc.a -ldl -lz -pthread > ${b}/${name}.build.log 2>&1 \
        || { echo "BUILD FAILED: ${name}"; tail -30 ${b}/${name}.build.log; } ) &
done < "${jd}/variants"
wait
mapfile -t V < <(awk 'NF { print $1 }' "${jd}/variants")
for v in "${V[@]}"; do [ -x ${b}/${v}.bin ] || exit 1; done
echo "built: ${V[*]}"

bind_of() # The numactl prefix of the memory policy `$1`.
{
    case "$1" in
        local)      echo "numactl --cpunodebind=${s0} --membind=${s0}" ;;
        interleave) echo "numactl --cpunodebind=${s0} --interleave=all" ;;
        all)        echo "numactl --interleave=all" ;;
        preferred)  echo "numactl --cpunodebind=${s0} --preferred=${s0%%,*}" ;;
        none)       echo "" ;;
    esac
}

echo "== runs"
last=0
for r in $(seq 1 ${ROUNDS}); do
    while read -r label policy threads args; do
        [ -z "${label}" ] && continue
        mkdir -p "${jd}/out/${label}"
        n=${#V[@]}
        for i in $(seq 0 $((n - 1))); do
            v=${V[$(( (i + r - 1) % n ))]}     # Rotate the variant order by round.
            dst="${jd}/out/${label}/${v}.r${r}.tsv"
            [ -s "${dst}" ] && continue
            if [ $(( $(date +%s) - t0 + last )) -gt $((BUDGET_MIN * 60)) ]; then
                echo "out of time before ${label}/${v}.r${r}; resubmit to continue"
                exit 0
            fi
            s=$(date +%s)
            $(bind_of ${policy}) ${b}/${v}.bin --threads ${threads} --reps ${REPS} ${args} > "${dst}.part" 2> "${dst}.err" < /dev/null \
                && mv "${dst}.part" "${dst}" || echo "FAILED ${label}/${v}.r${r}: $(tail -3 ${dst}.err)"
            last=$(( $(date +%s) - s ))
            echo "${label} ${v} r${r} (${last} s): build $(grep -P '^build\.mops\t' ${dst} 2> /dev/null | cut -f2 | paste -sd' ') q0 $(grep -P '^query0\.mops\t' ${dst} 2> /dev/null | cut -f2 | paste -sd' ') q1 $(grep -P '^query1\.mops\t' ${dst} 2> /dev/null | cut -f2 | paste -sd' ')"
        done
    done < "${jd}/configs"
done
echo "complete"
