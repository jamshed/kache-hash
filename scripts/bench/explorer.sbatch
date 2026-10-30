#!/bin/bash
#SBATCH --job-name=kh-bench
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --constraint=cascadelake
#SBATCH --exclude=d3230,d3231,d3204,d3233,d4041,d4100,d3091,d3092,d3093,d3094
#SBATCH --output=bench-%j.out

# Runs a plan of benchmarks on a whole node of explorer: builds every benchmark natively in the job, from a copy of the
# repository it is submitted from, runs `run.sh` for each line of the plan, and collates the outputs. A plan's line,
# which continues onto the next if it ends in `\`, is the variables of a `run.sh` matrix, as assignments, and then its
# workloads. Besides those of `run.sh`, the variables
# include `EXPERIMENT`, the line's name; `MEMORY`: by default, `LAUNCH` binds to the CPUs of socket 0, with memory
# preferred there—or with `MEMORY=interleave`, interleaved over all sockets, for tables beyond a socket's memory—unless
# `PLACEMENT` places each run instead (see `run.sh`); and `TIMEOUT`, a limit of each run, e.g. `3h`. The environment
# selects what of the plan runs, if set: `EXPERIMENTS`, the names of the lines to run; `TOOLS`, the tables to run (see
# `run.sh`); `WORKLOADS`, the workloads to run; and `SELECT`, the runs to run, by name (see `run.sh`). The job exits
# unless the node's CPU model contains `CPU_MODEL` (default: Platinum 8276).
# Usage, from the repository's root: sbatch scripts/bench/explorer.sbatch <data-dir> <out-dir> <plan>

data_dir="${1}"
out_dir="${2}"
plan="${3}"

if [ $# -lt 3 ]; then
    echo "Usage: sbatch $0 <data-dir> <out-dir> <plan>" >&2
    exit 1
fi

module load cmake/3.30.2
export LIBCLANG_PATH=/shared/EL9/explorer/anaconda3/2024.06/lib # CBL's bindings need a libclang; the system has none.
export BINDGEN_EXTRA_CLANG_ARGS="-isystem $(gcc -print-file-name=include)"  # That one lacks its builtin headers.

cpu_model="${CPU_MODEL:-Platinum 8276}"
model=$(lscpu | grep "^Model name")
[[ "${model}" == *"${cpu_model}"* ]] || { echo "Not a ${cpu_model} node: ${model}"; exit 3; }

build="/tmp/kache-hash-bench-${SLURM_JOB_ID}"   # Node-local, as concurrent jobs would clash in a shared one.
mkdir -p "${out_dir}" "${build}"
{ hostname; lscpu; numactl -H; cat /sys/kernel/mm/transparent_hugepage/enabled; } > "${out_dir}/node-${SLURM_JOB_ID}"

echo "Building in ${build}"
if ! rsync -a --exclude "build*/" --exclude external/ --exclude data/ --exclude local/ \
        "${SLURM_SUBMIT_DIR}/" "${build}/src/"; then
    echo "Failed copying ${SLURM_SUBMIT_DIR}" >&2
    exit 1
fi

if ! { cmake -S "${build}/src" -B "${build}/build" -DCMAKE_BUILD_TYPE=Release -DCBL=ON -DSSHASH=ON &&
        cmake --build "${build}/build" -j "$(nproc)"; } > "${out_dir}/build-${SLURM_JOB_ID}.log" 2>&1; then
    echo "Failed building; see ${out_dir}/build-${SLURM_JOB_ID}.log" >&2
    exit 1
fi

s0=$(lscpu -p=NODE,SOCKET | grep -v "^#" | awk -F, '$2 == 0 { print $1 }' | sort -un | paste -sd,)
status=0
while read -r line; do
    while [[ "${line}" == *\\ ]] && read -r next; do
        line="${line%\\} ${next}"
    done

    [[ "${line}" =~ ^[[:space:]]*(#|$) ]] && continue

    (
        eval "set -- ${line}"
        workloads=()
        for arg in "$@"; do
            if [[ "${arg}" =~ ^[A-Za-z_][A-Za-z0-9_]*= ]]; then
                export "${arg}"
            elif [ -z "${WORKLOADS:-}" ] || [[ " ${WORKLOADS} " == *" ${arg} "* ]]; then
                workloads+=("${arg}")
            fi
        done

        if { [ -n "${EXPERIMENTS:-}" ] && [[ " ${EXPERIMENTS} " != *" ${EXPERIMENT:-} "* ]]; } ||
                (( ${#workloads[@]} == 0 )); then
            echo "Skipping ${line}"
            exit 0
        fi

        echo "Running ${line}"
        if [ -n "${PLACEMENT:-}" ]; then
            export LAUNCH="${LAUNCH-${TIMEOUT:+timeout ${TIMEOUT}}}"
        else
            memory="--preferred=${s0%%,*}"
            [ "${MEMORY:-local}" == "interleave" ] && memory="--interleave=all"
            export LAUNCH="${LAUNCH-${TIMEOUT:+timeout ${TIMEOUT} }numactl --cpunodebind=${s0} ${memory}}"
        fi

        "${build}/src/scripts/bench/run.sh" "${build}/build/benchmark" "${data_dir}" "${out_dir}" "${workloads[@]}"
    ) < /dev/null || status=1
done < "${plan}"

python3 "${build}/src/scripts/bench/collate.py" "${out_dir}" || status=1

rm -rf "${build}"
exit ${status}
