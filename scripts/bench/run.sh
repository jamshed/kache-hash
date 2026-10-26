#!/bin/bash

# Runs the benchmark drivers built in a directory over workloads, a process per run, and writes each run's output to
# `<out-dir>/<workload>/<table>.<set|map>.<mt|seq>.<sizing>.T<threads>.r<rep>.tsv`. Runs with an output are skipped,
# so a rerun continues where an earlier one stopped. The environment sets the matrix:
#   TABLES:  the tables to run, of those built (default: kache-hash boost libcuckoo iceberg cbl sshash);
#   KINDS:   `set` and `map` (default: both);
#   SIZINGS: `presized`, for the input's distinct k-mers, and `half`, for half of them (default: both);
#   THREADS: the thread counts (default: 16); the sequential variants run at 1;
#   THREADS_<table>: the thread counts of one table, e.g. `THREADS_sshash`, instead of THREADS;
#   REPS:    the repetitions of each run (default: 3);
#   LAUNCH:  a prefix of each run, e.g. `numactl --cpunodebind=0 --membind=0` (default: none).
# The tables skip the configurations they do not support. Per input, the distinct k-mer count comes from a kache-hash
# set run grown from a small capacity, and the SSHash index from `sshash-index.sh`, both with the most threads of
# THREADS; both are kept in <data-dir>.

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
LAUNCH="${LAUNCH:-}"

max_t=$(printf "%s\n" ${THREADS} | sort -n | tail -n 1)
failed=0


# Sets `index` to the input that the workload $1 builds on, and `queries` to the inputs that it queries.
workload()
{
    case "${1}" in
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


# Sets `distinct` to the distinct k-mer count of the workload's input, counting it if not kept yet.
count_distinct()
{
    local kept="${data_dir}/${base}.distinct"

    if [ ! -s "${kept}" ]; then
        echo "Counting the distinct k-mers of ${index}"
        if ! ${LAUNCH} "${bin_dir}/bench_kache_hash" --index "${data_dir}/${index}" --set --threads "${max_t}" \
                --init-cap 1048576 > "${out}/distinct.tsv"; then
            echo "Failed counting the distinct k-mers of ${index}" >&2
            exit 1
        fi

        awk -F'\t' '$1 == "table.size" { print $2 }' "${out}/distinct.tsv" > "${kept}"
    fi

    distinct=$(cat "${kept}")
}


# Builds the SSHash index of the workload's input, if not kept yet, and copies its cost to the output.
index_sshash()
{
    if [ ! -s "${data_dir}/${base}.sshash" ]; then
        if ! GGCAT="${bin_dir}/ggcat" SSHASH="${bin_dir}/sshash" "$(dirname "${0}")/sshash-index.sh" \
                "${data_dir}/${index}" "${data_dir}/${base}" "${max_t}"; then
            echo "Failed building the SSHash index of ${index}" >&2
            exit 1
        fi
    fi

    if [ -s "${data_dir}/${base}.cost.tsv" ]; then
        cp "${data_dir}/${base}.cost.tsv" "${out}/sshash-build.tsv"
    fi
}


# Prints the run configurations of the matrix, as lines of table, kind, variant, sizing, and threads.
configs()
{
    local table kind threads t variant sizings sizing

    for table in ${tables}; do
        case "${table}" in
            cbl | sshash) sizings="na" ;;   # They can not be sized.
            iceberg)    # Its resizes duplicate and lose keys, even with one user thread.
                sizings=""
                [[ " ${SIZINGS} " == *" presized "* ]] && sizings="presized" ;;
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


# Runs the table $1 for the kind $2, the variant $3, the sizing $4, and the thread count $5, as the repetition $6.
run()
{
    local table="${1}" kind="${2}" variant="${3}" sizing="${4}" t="${5}" rep="${6}"
    local name="${table}.${kind}.${variant}.${sizing}.T${t}.r${rep}" args=() q

    [ -s "${out}/${name}.tsv" ] && return

    if [ "${table}" == "sshash" ]; then
        args+=(--index "${data_dir}/${base}.sshash")
    else
        args+=(--index "${data_dir}/${index}")
    fi

    for q in ${queries}; do
        args+=(--query "${data_dir}/${q}")
    done

    args+=(--threads "${t}")
    [ "${kind}" == "set" ] && args+=(--set)
    [ "${variant}" == "seq" ] && args+=(--seq)
    [ "${sizing}" == "presized" ] && args+=(--init-cap "${distinct}")
    [ "${sizing}" == "half" ] && args+=(--init-cap "$(( distinct / 2 ))")

    if ${LAUNCH} "${bin_dir}/bench_${table//-/_}" "${args[@]}" > "${out}/${name}.tsv.tmp" 2> "${out}/${name}.err"; then
        mv "${out}/${name}.tsv.tmp" "${out}/${name}.tsv"
        [ -s "${out}/${name}.err" ] || rm -f "${out}/${name}.err"
        echo "Ran ${wl}/${name}"
    elif grep -q "^Unsupported configuration" "${out}/${name}.err"; then
        rm -f "${out}/${name}.tsv.tmp" "${out}/${name}.err"
    else
        echo "Failed running ${wl}/${name}; see ${out}/${name}.err" >&2
        failed=$(( failed + 1 ))
    fi
}


tables=""
for table in ${TABLES}; do
    if [ -x "${bin_dir}/bench_${table//-/_}" ]; then
        tables="${tables} ${table}"
    else
        echo "Skipping ${table}: no driver in ${bin_dir}"
    fi
done

for wl in "$@"; do
    workload "${wl}"
    base="${index%.gz}"
    base="${base%.fa}"
    base="${base%.fastq}"
    out="${out_dir}/${wl}"
    mkdir -p "${out}"

    if configs | grep -qE " (presized|half) "; then
        count_distinct
    fi

    if configs | grep -q "^sshash "; then
        index_sshash
    fi

    for rep in $(seq 1 "${REPS}"); do   # Repetitions interleave, to spread drifts over configurations.
        while read -r table kind variant sizing t; do
            run "${table}" "${kind}" "${variant}" "${sizing}" "${t}" "${rep}" < /dev/null
        done < <(configs)
    done
done

if (( failed > 0 )); then
    echo "${failed} runs failed" >&2
    exit 1
fi
