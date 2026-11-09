#!/usr/bin/env bash
# Builds test/check.cpp, test/describers.cpp, test/counter.cpp with the counting library, test/oriented.cpp with the
# batch and the lane describers, and test/capi.cpp with the C interface library, locally, as x86_64 (under Rosetta on
# Apple silicon), and checks that the C interface's header is C. Runs the describers' check on E. coli, chr21 and 1M
# SRR5833294 reads, the counting library's checks, the oriented operations' and the C interface's checks on E. coli
# and chr21 at all threads, then the tables' correctness matrix: E. coli and chr21, each against the other as
# negatives; sets, and maps of 64-, 32- and 8-bit counts; sequential, and concurrent at 1 and all threads; presized,
# and resizing from 100,000 and from 10 k-mers; whole fragments, raw sequences, and ragged batches; and concurrent at
# all threads with bucket counts rounded up to powers of 2, presized and from 10 k-mers, over whole fragments.
# Usage: test/check.sh [quick] [extra compile flags...]; quick skips chr21 and the reads.
set -u
root="$(cd "$(dirname "$0")/.." && pwd)"
quick=0; [ "${1:-}" = quick ] && { quick=1; shift; }
out="${root}/build-local/check$(echo "$*" | tr -c 'A-Za-z0-9=' '_')"; mkdir -p "${out}"   # Per flags, so concurrent runs don't share a binary.
for t in check describers counter oriented oriented-lanes capi; do
    srcs="${root}/test/${t%-lanes}.cpp" defs=""
    [ ${t} = counter ] && srcs+=" ${root}/counter/Fastx_Reader.cpp ${root}/counter/Kmer_Set_Counter.cpp"
    [ ${t} = oriented-lanes ] && defs="-DKACHE_HASH_LANES=1"
    [ ${t} = capi ] && srcs+=" ${root}/capi/*.cpp"
    clang++ -arch x86_64 -std=c++20 -O2 -march=haswell -Wall -Wextra -Wno-unknown-attributes -I "${root}/include" \
        -I "${root}/support" -I "${root}/counter" -I "${root}/capi" -pthread ${defs} "$@" ${srcs} -o "${out}/${t}" -lz || exit 1
done
clang -arch x86_64 -std=c11 -Wall -Wextra -Wpedantic -Werror -fsyntax-only -x c "${root}/capi/kache_capi.h" || exit 1
T=$(sysctl -n hw.ncpu)
fail=0
inputs="${root}/data/ecoli.fa.gz"; [ ${quick} = 0 ] && inputs+=" ${root}/data/chr21.fa.gz ${root}/data/SRR5833294.1M.fastq.gz"
"${out}/describers" ${inputs} || fail=$((fail + 1))
"${out}/counter" check "${out}/counter-data" ${T} || fail=$((fail + 1))
for g in ecoli chr21; do
    neg=chr21; [ ${g} = chr21 ] && neg=ecoli
    [ ${quick} = 1 ] && [ ${g} = chr21 ] && continue
    for t in oriented oriented-lanes capi; do
        "${out}/${t}" "${root}/data/${g}.fa.gz" "${root}/data/${neg}.fa.gz" ${T} || fail=$((fail + 1))
    done
done
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
    for cap in 0 10; do
        for kind in set map map32 map8; do
            run "${root}/data/${g}.fa.gz" "${root}/data/${neg}.fa.gz" mt-${kind} ${T} ${cap} stream pow2
        done
    done
done
echo "failures: ${fail}"
exit ${fail}
