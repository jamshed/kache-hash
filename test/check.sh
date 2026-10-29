#!/usr/bin/env bash
# Builds test/check.cpp and test/describers.cpp locally, as x86_64 (under Rosetta on Apple silicon), and runs the
# describers' check on E. coli, chr21 and 1M SRR5833294 reads, then the tables' correctness matrix: E. coli and chr21,
# each against the other as negatives; sets, and maps of 64-, 32- and 8-bit counts; sequential, and concurrent at 1 and
# all threads; presized, and resizing from 100,000 and from 10 k-mers; whole fragments, raw sequences, and ragged
# batches.
# Usage: test/check.sh [quick] [extra compile flags...]; quick skips chr21 and the reads.
set -u
root="$(cd "$(dirname "$0")/.." && pwd)"
quick=0; [ "${1:-}" = quick ] && { quick=1; shift; }
out="${root}/build-local/check$(echo "$*" | tr -c 'A-Za-z0-9=' '_')"; mkdir -p "${out}"   # Per flags, so concurrent runs don't share a binary.
for t in check describers; do
    clang++ -arch x86_64 -std=c++20 -O2 -march=haswell -DDEVELOP -DL1_CACHE_LINE_SIZE=64 -DXXH_INLINE_ALL -DXXH_STATIC_LINKING_ONLY \
        -DXXH_FORCE_MEMORY_ACCESS=0 -Wall -Wextra -Wno-unknown-attributes -I "${root}/include" -pthread "$@" \
        "${root}/test/${t}.cpp" -o "${out}/${t}" -lz || exit 1
done
T=$(sysctl -n hw.ncpu)
fail=0
inputs="${root}/data/ecoli.fa.gz"; [ ${quick} = 0 ] && inputs+=" ${root}/data/chr21.fa.gz ${root}/data/SRR5833294.1M.fastq.gz"
"${out}/describers" ${inputs} || fail=$((fail + 1))
run() { "${out}/check" "$@" || fail=$((fail + 1)); }
for g in ecoli chr21; do
    neg=chr21; [ ${g} = chr21 ] && neg=ecoli
    [ ${quick} = 1 ] && [ ${g} = chr21 ] && continue
    for cap in 0 100000 10; do
        for mode in stream raw multi; do
            for kind in set map map32 map8; do
                run "${root}/data/${g}.fa.gz" "${root}/data/${neg}.fa.gz" seq-${kind} 1 ${cap} ${mode}
                for t in 1 ${T}; do
                    run "${root}/data/${g}.fa.gz" "${root}/data/${neg}.fa.gz" mt-${kind} ${t} ${cap} ${mode}
                done
            done
        done
    done
done
echo "failures: ${fail}"
exit ${fail}
