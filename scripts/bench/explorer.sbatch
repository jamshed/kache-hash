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

# Runs `run.sh` on a whole Cascade Lake node of explorer: builds every benchmark natively in the job, from a copy of the
# repository it is submitted from, runs them, and collates the outputs. The matrix variables of `run.sh` pass through;
# `LAUNCH` defaults to the CPUs of socket 0, with memory preferred there.
# Usage, from the repository's root: sbatch scripts/bench/explorer.sbatch <data-dir> <out-dir> <workload>...

data_dir="${1}"
out_dir="${2}"

if [ $# -lt 3 ]; then
    echo "Usage: sbatch $0 <data-dir> <out-dir> <workload>..." >&2
    exit 1
fi

shift 2

module load cmake/3.30.2
export LIBCLANG_PATH=/shared/EL9/explorer/anaconda3/2024.06/lib # CBL's bindings need a libclang; the system has none.
export BINDGEN_EXTRA_CLANG_ARGS="-isystem $(gcc -print-file-name=include)"  # That one lacks its builtin headers.

lscpu | grep -q "Platinum 8276" || { echo "Not a Cascade Lake node: $(lscpu | grep '^Model name')"; exit 3; }

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
export LAUNCH="${LAUNCH-numactl --cpunodebind=${s0} --preferred=${s0%%,*}}"

"${build}/src/scripts/bench/run.sh" "${build}/build/benchmark" "${data_dir}" "${out_dir}" "$@"
status=$?
python3 "${build}/src/scripts/bench/collate.py" "${out_dir}" || status=1

rm -rf "${build}"
exit ${status}
