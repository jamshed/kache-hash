#!/bin/bash

# Builds the SSHash index of the k-mers of an input, as SSHash requires: GGCAT's eulertigs of the input, keeping every
# k-mer, and then the SSHash dictionary over them. Writes `<out-prefix>.eulertigs.fa`, `<out-prefix>.sshash`, and the
# cost of each stage—wall time, CPU time, and max RSS—to `<out-prefix>.cost.tsv`. The minimizer length m defaults to
# SSHash's lower bound for its input. The tools are `GGCAT` and `SSHASH` if set, or else `ggcat` and `sshash` from
# `PATH`.

input="${1}"
out_prefix="${2}"
thread_c="${3}"
m="${4}"

if [ $# -lt 3 ]; then
    echo "Usage: $0 <input> <out-prefix> <threads> <(optional) minimizer-length>" >&2
    exit 1
fi

ggcat="${GGCAT:-ggcat}"
sshash="${SSHASH:-sshash}"
k=31
tigs="${out_prefix}.eulertigs.fa"
index="${out_prefix}.sshash"
cost="${out_prefix}.cost.tsv"
tmp_dir="${out_prefix}.tmp"


# Runs the command $2... as the stage $1 under GNU `time`, and appends its cost to the cost file.
run_stage()
{
    local stage="${1}" log="${tmp_dir}/${1}.time"
    shift

    echo "Running ${stage}"
    if ! /usr/bin/time -v -o "${log}" "$@" > "${tmp_dir}/${stage}.log" 2>&1; then
        tail "${tmp_dir}/${stage}.log" >&2
        echo "Failed running ${stage}" >&2
        exit 1
    fi

    awk -F': ' -v s="${stage}" '
        /Elapsed \(wall clock\)/ {    # As [h:]m:s.
            n = split($2, t, ":"); w = 0
            for(i = 1; i <= n; ++i) w = w * 60 + t[i]
            print s ".seconds\t" w }
        /User time/ { cpu += $2 }
        /System time/ { cpu += $2 }
        /Maximum resident set size/ { print s ".max_rss_bytes\t" $2 * 1024 }
        END { print s ".cpu_seconds\t" cpu }' "${log}" >> "${cost}"
}


rm -rf "${tmp_dir}"
mkdir -p "${tmp_dir}"
printf "input\t%s\nk\t%s\nthreads\t%s\n" "${input}" "${k}" "${thread_c}" > "${cost}"

# The minimum multiplicity 1 keeps every k-mer; eulertigs hold each one once.
run_stage ggcat "${ggcat}" build -k "${k}" -s 1 --eulertigs -j "${thread_c}" -t "${tmp_dir}/ggcat" \
    -o "${tigs}" "${input}"
rm -f "${out_prefix}.eulertigs.stats.log"   # GGCAT's per-second progress counters.

bases=$(grep -v '^>' "${tigs}" | tr -d '\n' | wc -c | tr -d ' ')
if [ -z "${m}" ]; then  # SSHash's lower bound on m: ceil(log_4(N)) + 1, for N input bases.
    e=0
    while (( 4 ** e < bases )); do (( e++ )); done
    m=$(( e + 1 ))
fi

printf "eulertigs.bases\t%s\nsshash.m\t%s\n" "${bases}" "${m}" >> "${cost}"

run_stage sshash "${sshash}" build -i "${tigs}" -k "${k}" -m "${m}" -t "${thread_c}" -d "${tmp_dir}/sshash" \
    -o "${index}"
printf "sshash.index_bytes\t%s\n" "$(wc -c < "${index}" | tr -d ' ')" >> "${cost}"

rm -rf "${tmp_dir}"
