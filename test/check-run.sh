#!/bin/bash
# Builds test/check.cpp, test/describers.cpp, test/oriented.cpp with the node's describer and with the batch one, and
# test/capi.cpp with the C interface library, natively from a source tree pushed to explorer (scripts/ab/push.sh),
# with extra compile flags, and runs the describers' check on E. coli, chr21 and 1M SRR5833294 reads, the oriented
# operations' and the C interface's checks on E. coli at `T` threads, then check.sh's matrix for E. coli on the node,
# with chr21 as negatives, at 1 and `T` threads, for sets and maps of 64-, 32- and 8-bit counts, with its cases of
# bucket counts rounded up to powers of 2 at `T` threads.
# Usage: check-run.sh <source> <T> [extra compile flags...], e.g. as an A/B job's KH_PRE. KH_DATA overrides the data
# directory, e.g. for copies staged in node-local /tmp; KH_NO_DESC=1 skips the describers' check.
set -u
src="$1"; T="$2"; shift 2
data=${KH_DATA:-/scratch/md.khan/kache-hash/data}
b=/tmp/khc-${SLURM_JOB_ID:-$$}/$(echo "x $*" | tr -c 'A-Za-z0-9=' '_')
mkdir -p ${b}
rsync -a --exclude "build*/" --exclude external/ ~/kh-perf/src/${src}/ ${b}/src/
for t in check describers oriented oriented-batch capi; do
    srcs="${b}/src/test/${t%-batch}.cpp" defs=""
    [ ${t} = oriented-batch ] && defs="-DKACHE_HASH_LANES=0"
    [ ${t} = capi ] && srcs+=" ${b}/src/capi/*.cpp"
    g++ -std=c++20 -O2 -march=native -DDEVELOP -DL1_CACHE_LINE_SIZE=64 -DXXH_INLINE_ALL -DXXH_STATIC_LINKING_ONLY \
        -DXXH_FORCE_MEMORY_ACCESS=0 -Wall -Wextra ${defs} "$@" -I ${b}/src/include -I ${b}/src/capi -pthread ${srcs} -o ${b}/${t} -lz \
        || { echo "check [${src} $*]: BUILD FAILED"; exit 1; }
done
fail=0; c=0
[ "${KH_NO_DESC:-0}" = 1 ] || ${b}/describers ${data}/ecoli.fa.gz ${data}/chr21.fa.gz ${data}/SRR5833294.1M.fastq.gz || fail=$((fail + 1))
for t in oriented oriented-batch capi; do
    c=$((c + 1)); ${b}/${t} ${data}/ecoli.fa.gz ${data}/chr21.fa.gz ${T} || fail=$((fail + 1))
done
run() { c=$((c + 1)); ${b}/check "$@" > ${b}/last.txt 2>&1 || { fail=$((fail + 1)); cat ${b}/last.txt; }; }
for cap in 0 100000 10; do
    for mode in stream raw multi; do
        for kind in set map map32 map8; do
            run ${data}/ecoli.fa.gz ${data}/chr21.fa.gz seq-${kind} 1 ${cap} ${mode}
            for t in 1 ${T}; do
                run ${data}/ecoli.fa.gz ${data}/chr21.fa.gz mt-${kind} ${t} ${cap} ${mode}
            done
        done
    done
done
for cap in 0 10; do
    for kind in set map map32 map8; do
        run ${data}/ecoli.fa.gz ${data}/chr21.fa.gz mt-${kind} ${T} ${cap} stream pow2
    done
done
echo "check [${src} $*] T=${T}: ${c} cases, failures: ${fail}"
exit ${fail}
