#!/bin/bash

# Runs the benchmark drivers built in a directory over workloads, a process per run, and writes each run's output to
# `<out-dir>/<workload>/<table>.<set|map>.<mt|seq>.<sizing>.T<threads>.r<rep>.tsv`. Runs with an output are skipped,
# so a rerun continues where an earlier one stopped. The environment sets the matrix:
#   TABLES:  the tables to run, of those built (default: kache-hash boost libcuckoo iceberg cbl sshash);
#   TOOLS:   if set, the only tables of TABLES to run, to select tables across a plan's matrices;
#   KINDS:   `set` and `map` (default: both);
#   SIZINGS: `presized`, for the input's distinct k-mers, and `half`, for half of them (default: both);
#   THREADS: the thread counts (default: 16); the sequential variants run at 1;
#   THREADS_<table>: the thread counts of one table, e.g. `THREADS_sshash`, instead of THREADS;
#   QUERY_THREADS: further thread counts to query each table with, those below its run's, reported as the phases
#            `query<i>.T<count>` (default: none);
#   QUERIES: if set, the inputs to query instead of the workload's, e.g. empty to only build;
#   REPS:    the repetitions of each run (default: 3);
#   SELECT:  if set, an extended regular expression that a run's name—its output's path in `<out-dir>`, less `.tsv`—must
#            match for it to run, to select runs across a plan's matrices, e.g. `^reference/kache-hash\.` or `\.T56\.`;
#   PLACEMENT: if set, `spread` or `compact`: each run takes a core per thread, alternating between the sockets or
#            filling one socket before the next, with its memory interleaved over the NUMA nodes (default: none);
#   LAUNCH:  a prefix of each run, ahead of PLACEMENT's, e.g. `timeout 3h`, or without PLACEMENT, a binding such as
#            `numactl --cpunodebind=0 --membind=0` (default: none);
#   SSHASH_THREADS: the thread count of SSHash's index build (default: the most of THREADS);
#   SSHASH_M: the minimizer length of SSHash (default: `sshash-index.sh`'s);
#   VERSION: the code's version, which each run's output records as `version` (default: the repository's version, as
#            `git describe --always --dirty` names it, with a hash of any uncommitted changes; or `unknown`).
# The tables skip the configurations they do not support. Per workload, as its runs need them, the distinct k-mer count
# comes from a kache-hash set run grown from a small capacity, with the most threads of THREADS, and the SSHash index
# from `sshash-index.sh`; both launch as the runs do, and are kept with the outputs, as `distinct.tsv` and `index.*`, so
# that each machine measures its own. Each launch's command line is kept beside its output, as `<run>.launch`.

bin_dir="${1}"
data_dir="${2}"
out_dir="${3}"

if [ $# -lt 4 ]; then
    echo "Usage: $0 <bin-dir> <data-dir> <out-dir> <workload>..." >&2
    exit 1
fi

shift 3

TABLES="${TABLES:-kache-hash boost libcuckoo iceberg cbl sshash}"
KINDS="${KINDS:-set map}"
SIZINGS="${SIZINGS:-presized half}"
THREADS="${THREADS:-16}"
REPS="${REPS:-3}"
SELECT="${SELECT:-}"
PLACEMENT="${PLACEMENT:-}"
LAUNCH="${LAUNCH:-}"

max_t=$(printf "%s\n" ${THREADS} | sort -n | tail -n 1)
SSHASH_THREADS="${SSHASH_THREADS:-${max_t}}"
failed=0


# Prints the CPUs that runs take, a core each, in the order of PLACEMENT: `spread` alternates between the sockets, and
# `compact` takes a socket's cores before the next one's.
placement_cpus()
{
    local range cpu topo

    for range in $(awk '/^Cpus_allowed_list/ { gsub(",", " ", $2); print $2 }' /proc/self/status); do
        for cpu in $(seq "${range%-*}" "${range#*-}"); do
            topo="/sys/devices/system/cpu/cpu${cpu}/topology"
            echo "$(cat "${topo}/physical_package_id") $(cat "${topo}/core_id") ${cpu}"
        done
    done | sort -k1,1n -k3,3n | awk -v placement="${PLACEMENT}" '
        ($1, $2) in taken { next }  # Another hardware thread of a core taken.
        { taken[$1, $2]; if (!($1 in n)) sockets[s++] = $1; cpus[$1, n[$1]++] = $3 }
        END {
            for (j = 0; j < s; ++j) if (n[sockets[j]] > m) m = n[sockets[j]]
            if (placement == "compact") {
                for (j = 0; j < s; ++j) for (i = 0; i < n[sockets[j]]; ++i) print cpus[sockets[j], i]
            } else {
                for (i = 0; i < m; ++i) for (j = 0; j < s; ++j) if (i < n[sockets[j]]) print cpus[sockets[j], i]
            }
        }'
}


# Sets `placement` to the prefix that places a launch of $1 threads as PLACEMENT sets, if set: on the first $1 of its
# CPUs, with memory interleaved over all the NUMA nodes.
place()
{
    placement=""
    if [ -n "${PLACEMENT}" ]; then
        placement="numactl --interleave=all --physcpubind=$(IFS=,; echo "${cpus[*]:0:${1}}")"
    fi
}


# Runs the command $3... as a launch of $2 threads, prefixed with LAUNCH and placed as PLACEMENT sets, and writes its
# command line into the file $1.
launch()
{
    local file="${1}" t="${2}"
    shift 2

    place "${t}"
    echo ${LAUNCH} ${placement} "$@" > "${file}"
    ${LAUNCH} ${placement} "$@"
}


# Sets `index` to the input that the workload $1 builds on, and `queries` to the inputs that it queries.
workload()
{
    case "${1}" in
        # The paper's workloads: each index queried with k-mers it holds (positive), k-mers of another domain of life
        # (negative), and an application's input (mixed).
        reference)
            index=CHM13.fa.gz
            queries="CHM13.fa.gz SRR5901135.fastq.gz SRR37356339.4M.fastq.gz" ;;
        short-reads)
            index=SRR37356339.20M.fastq.gz
            queries="SRR37356339.4M.fastq.gz SRR5901135.fastq.gz CHM13.fa.gz" ;;
        long-reads)
            index=SRR10382246.360K.fastq.gz
            queries="SRR10382246.45K.fastq.gz SRR5901135.fastq.gz SRR37356339.4M.fastq.gz" ;;
        pangenome)
            index=ecoli-3682.fa.gz
            queries="ecoli-3682.every8.fa.gz SRR37356339.4M.fastq.gz SRR5901135.fastq.gz" ;;
        metagenome)
            index=humgut-1in20.fa.gz
            queries="humgut-1in20.every8.fa.gz SRR37356339.4M.fastq.gz ERR321482.fastq.gz" ;;
        # The earlier campaigns' workloads.
        grch38)
            index=GRCh38.fa.gz
            queries="GRCh38.fa.gz ecoli.fa.gz SRR5833294.fastq.gz SRR5901135_1.fastq.gz SRR5901135_2.fastq.gz" ;;
        srr5833294)
            index=SRR5833294.fastq.gz
            queries="SRR5833294.fastq.gz SRR5901135_1.fastq.gz" ;;
        chr1)
            index=chr1.fa.gz
            queries="chr1.fa.gz ecoli.fa.gz" ;;
        srr5901135)
            index=SRR5901135_1.fastq.gz
            queries="SRR5901135_1.fastq.gz ecoli.fa.gz" ;;
        ecoli-1024)
            index=ecoli-1024.fa.gz
            queries="ecoli.fa.gz SRR5901135_1.fastq.gz chr21.fa.gz" ;;
        # The small workloads, for local runs.
        chr21)
            index=chr21.fa.gz
            queries="chr21.fa.gz ecoli.fa.gz" ;;
        srr5833294-1m)
            index=SRR5833294.1M.fastq.gz
            queries="SRR5833294.1M.fastq.gz SRR5901135_1.500K.fastq.gz" ;;
        ecoli-32)
            index=ecoli-32.fa.gz
            queries="ecoli.fa.gz SRR5901135_1.500K.fastq.gz chr21.fa.gz" ;;
        *)
            echo "Unknown workload ${1}" >&2
            exit 1 ;;
    esac
}


# Sets `distinct` to the distinct k-mer count of the workload's input, counting it if not counted yet.
count_distinct()
{
    if [ ! -s "${out}/distinct.tsv" ]; then
        echo "Counting the distinct k-mers of ${index}"
        if ! launch "${out}/distinct.launch" "${max_t}" "${bin_dir}/bench_kache_hash" --index "${data_dir}/${index}" \
                --set --threads "${max_t}" --init-cap 1048576 > "${out}/distinct.tsv.tmp"; then
            echo "Failed counting the distinct k-mers of ${index}" >&2
            exit 1
        fi

        mv "${out}/distinct.tsv.tmp" "${out}/distinct.tsv"
    fi

    distinct=$(awk -F'\t' '$1 == "table.size" { print $2 }' "${out}/distinct.tsv")
}


# Builds the SSHash index of the workload's input, if not built yet: its cost's last line marks it complete. Sets
# `indexed` to whether it is built.
index_sshash()
{
    indexed=1
    if ! grep -q "^sshash.index_bytes" "${out}/index.cost.tsv" 2> /dev/null; then
        if ! launch "${out}/index.launch" "${SSHASH_THREADS}" env GGCAT="${bin_dir}/ggcat" \
                SSHASH="${bin_dir}/sshash" "$(dirname "${0}")/sshash-index.sh" "${data_dir}/${index}" "${out}/index" \
                "${SSHASH_THREADS}" "${SSHASH_M:-}"; then
            echo "Failed building the SSHash index of ${index}; skipping its SSHash runs" >&2
            failed=$(( failed + 1 ))
            indexed=0
        fi
    fi
}


# Prints the run configurations of the matrix, as lines of table, kind, variant, sizing, and threads.
configs()
{
    local table kind threads t variant sizings sizing

    for table in ${tables}; do
        case "${table}" in
            cbl | sshash) sizings="na" ;;   # They can not be sized.
            *) sizings="${SIZINGS}" ;;
        esac

        threads="THREADS_${table//-/_}"
        for kind in ${KINDS}; do
            for t in ${!threads:-${THREADS}}; do
                for variant in mt $([ "${t}" -eq 1 ] && echo seq); do
                    for sizing in ${sizings}; do
                        echo "${table} ${kind} ${variant} ${sizing} ${t}"
                    done
                done
            done
        done
    done
}


# Prints the runs of the workload `wl` that SELECT selects, as lines of name, table, kind, variant, sizing, and threads.
# Repetitions interleave, to spread drifts over configurations.
runs()
{
    local rep table kind variant sizing t name

    for rep in $(seq 1 "${REPS}"); do
        while read -r table kind variant sizing t; do
            name="${table}.${kind}.${variant}.${sizing}.T${t}.r${rep}"
            if [ -z "${SELECT}" ] || [[ "${wl}/${name}" =~ ${SELECT} ]]; then
                echo "${name} ${table} ${kind} ${variant} ${sizing} ${t}"
            fi
        done < <(configs)
    done
}


# Runs the table $2 for the kind $3, the variant $4, the sizing $5, and the thread count $6, as the run named $1.
run()
{
    local name="${1}" table="${2}" kind="${3}" variant="${4}" sizing="${5}" t="${6}" args=() q

    [ -s "${out}/${name}.tsv" ] && return

    if [ "${table}" == "sshash" ]; then
        args+=(--index "${out}/index.sshash")
    else
        args+=(--index "${data_dir}/${index}")
    fi

    for q in ${queries}; do
        args+=(--query "${data_dir}/${q}")
    done

    args+=(--threads "${t}")
    for q in ${QUERY_THREADS:-}; do
        [ "${q}" -lt "${t}" ] && args+=(--query-threads "${q}")
    done

    [ "${kind}" == "set" ] && args+=(--set)
    [ "${variant}" == "seq" ] && args+=(--seq)
    [ "${sizing}" == "presized" ] && args+=(--init-cap "${distinct}")
    [ "${sizing}" == "half" ] && args+=(--init-cap "$(( distinct / 2 ))")

    if launch "${out}/${name}.launch" "${t}" "${bin_dir}/bench_${table//-/_}" "${args[@]}" \
            > "${out}/${name}.tsv.tmp" 2> "${out}/${name}.err"; then
        printf "version\t%s\n" "${VERSION}" >> "${out}/${name}.tsv.tmp"
        mv "${out}/${name}.tsv.tmp" "${out}/${name}.tsv"
        [ -s "${out}/${name}.err" ] || rm -f "${out}/${name}.err"
        echo "Ran ${wl}/${name}"
    elif grep -q "^Unsupported configuration" "${out}/${name}.err"; then
        rm -f "${out}/${name}.tsv.tmp" "${out}/${name}.err" "${out}/${name}.launch"
    else
        echo "Failed running ${wl}/${name}; see ${out}/${name}.err" >&2
        failed=$(( failed + 1 ))
    fi
}


# Prints the version of the code, as `git describe` names that of this script's repository, with a hash of its
# uncommitted changes, if any; or `unknown`, outside a repository.
version()
{
    local repo describe hash

    repo="$(dirname "${0}")"
    if ! describe=$(git -C "${repo}" describe --always --dirty 2> /dev/null); then
        echo unknown
    elif [[ "${describe}" == *-dirty ]]; then
        hash=$(git -C "${repo}" diff HEAD --binary | { sha256sum 2> /dev/null || shasum -a 256; })
        echo "${describe}-${hash:0:8}"
    else
        echo "${describe}"
    fi
}


tables=""
for table in ${TABLES}; do
    if [ -n "${TOOLS:-}" ] && [[ " ${TOOLS} " != *" ${table} "* ]]; then
        continue
    elif [ -x "${bin_dir}/bench_${table//-/_}" ]; then
        tables="${tables} ${table}"
    else
        echo "Skipping ${table}: no driver in ${bin_dir}"
    fi
done

if [ -z "${tables}" ]; then
    echo "No tables to run"
    exit 0
fi

if [ -n "${SELECT}" ]; then
    [[ "" =~ ${SELECT} ]]
    if (( $? == 2 )); then
        echo "SELECT is not an extended regular expression: ${SELECT}" >&2
        exit 1
    fi
fi

if [ -n "${PLACEMENT}" ]; then
    if [[ " spread compact " != *" ${PLACEMENT} "* ]] || ! command -v numactl > /dev/null; then
        echo "PLACEMENT is spread or compact, and needs numactl" >&2
        exit 1
    fi

    cpus=($(placement_cpus))
    most=$(printf "%s\n" "${max_t}" "${SSHASH_THREADS}" $(configs | cut -d' ' -f5) | sort -n | tail -n 1)
    if (( most > ${#cpus[@]} )); then
        echo "Placing ${most} threads needs more than the ${#cpus[@]} cores available" >&2
        exit 1
    fi
fi

VERSION="${VERSION:-$(version)}"

for wl in "$@"; do
    workload "${wl}"
    queries="${QUERIES-${queries}}"
    out="${out_dir}/${wl}"

    selected=$(runs)
    if [ -z "${selected}" ]; then
        echo "Skipping ${wl}: no runs selected"
        continue
    fi

    mkdir -p "${out}"
    if grep -qE " (presized|half) " <<< "${selected}"; then
        count_distinct
    fi

    indexed=0
    if grep -q "^sshash\." <<< "${selected}"; then
        index_sshash
    fi

    while read -r name table kind variant sizing t; do
        [ "${table}" == "sshash" ] && (( ! indexed )) && continue
        run "${name}" "${table}" "${kind}" "${variant}" "${sizing}" "${t}" < /dev/null
    done <<< "${selected}"
done

if (( failed > 0 )); then
    echo "${failed} runs failed" >&2
    exit 1
fi
