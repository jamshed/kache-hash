#!/bin/bash
# Builds test/check.cpp and test/describers.cpp natively from a source tree pushed to explorer (scripts/ab/push.sh),
# with extra compile flags, and runs the describers' check on E. coli, chr21 and 1M SRR5833294 reads, then check.sh's
# matrix for E. coli on the node, with chr21 as negatives, at 1 and `T` threads.
# Usage: check-run.sh <source> <T> [extra compile flags...], e.g. as an A/B job's KH_PRE.
set -u
src="$1"; T="$2"; shift 2
data=/scratch/md.khan/kache-hash/data
b=/tmp/khc-${SLURM_JOB_ID:-$$}/$(echo "x $*" | tr -c 'A-Za-z0-9=' '_')
mkdir -p ${b}
rsync -a --exclude "build*/" --exclude external/ ~/kh-perf/src/${src}/ ${b}/src/
for t in check describers; do
    g++ -std=c++20 -O2 -march=native -DDEVELOP -DL1_CACHE_LINE_SIZE=64 -DXXH_INLINE_ALL -DXXH_STATIC_LINKING_ONLY \
        -DXXH_FORCE_MEMORY_ACCESS=0 -Wall -Wextra "$@" -I ${b}/src/include -pthread ${b}/src/test/${t}.cpp -o ${b}/${t} -lz \
        || { echo "check [${src} $*]: BUILD FAILED"; exit 1; }
done
fail=0; c=0
${b}/describers ${data}/ecoli.fa.gz ${data}/chr21.fa.gz ${data}/SRR5833294.1M.fastq.gz || fail=$((fail + 1))
run() { c=$((c + 1)); ${b}/check "$@" > ${b}/last.txt 2>&1 || { fail=$((fail + 1)); cat ${b}/last.txt; }; }
for cap in 0 100000 1000; do
    for mode in stream raw multi; do
        run ${data}/ecoli.fa.gz ${data}/chr21.fa.gz seq-set 1 ${cap} ${mode}
        run ${data}/ecoli.fa.gz ${data}/chr21.fa.gz seq-map 1 ${cap} ${mode}
        for t in 1 ${T}; do
            run ${data}/ecoli.fa.gz ${data}/chr21.fa.gz mt-set ${t} ${cap} ${mode}
            run ${data}/ecoli.fa.gz ${data}/chr21.fa.gz mt-map ${t} ${cap} ${mode}
        done
    done
done
echo "check [${src} $*] T=${T}: ${c} cases, failures: ${fail}"
exit ${fail}
