#!/bin/bash

# Feeds experiments of a plan to a node of explorer, a job at a time, each experiment's runs into `<out-root>/<name>`,
# as their runs' names do not tell the experiments apart: an experiment's job runs again, continuing where the last
# stopped (`run.sh` skips finished runs), until one completes, and then the next experiment's starts. The jobs build and
# run from the repository that this runs in, a snapshot keeping its `.git` (see `explorer.sbatch`). They run in
# `sharing`, for at most `TIME` (default: 00:59:00); or with `RESERVATION` set, in that reservation, in the partition
# `PARTITION` (default: reservation), for at most `TIME` (default: 1-00:00:00). A job that times out runs again if it
# wrote any output; one that finds the node busy (status 4) runs again 10 minutes later, up to 6 times; any other end
# stops the feeding, which a rerun resumes. The jobs get the environment's selection (`TOOLS`, `WORKLOADS`, `SELECT`)
# and settings (e.g. `TIMEOUT`, `STRAY_CPU`), but not the variables of `run.sh`'s matrices, which the plan's lines set
# alone. The feeding is logged to `<out-root>/feed-<node>.log`, and the experiments done listed in
# `<out-root>/feed-<node>.done`, which a rerun skips. Run it on a login node, e.g. under `nohup`, one per node, each
# with an out-root of its own.
# The partitions are explorer's, `sharing` and `reservation`: on another cluster, `where` below takes its own.
# Usage, from the repository's root: scripts/bench/feed.sh <node> <data-dir> <out-root> <plan> <experiment>...

node="${1}"
data_dir="${2}"
out_root="${3}"
plan="${4}"

if [ $# -lt 5 ]; then
    echo "Usage: $0 <node> <data-dir> <out-root> <plan> <experiment>..." >&2
    exit 1
fi

if [ ! -f scripts/bench/explorer.sbatch ]; then
    echo "Run $0 from the repository's root" >&2
    exit 1
fi

shift 4
busy_wait=600   # Seconds to wait for a busy node before trying it again.
busy_tries=6    # Most retries of a busy node.

if [ -n "${RESERVATION:-}" ]; then
    where=(--partition="${PARTITION:-reservation}" --reservation="${RESERVATION}" --time="${TIME:-1-00:00:00}")
else
    where=(--partition=sharing --time="${TIME:-00:59:00}")
fi

# The variables of `run.sh`'s matrices and of the plan's lines, unset for the jobs.
matrix=(EXPERIMENT LAUNCH MEMORY TABLES KINDS VARIANTS SIZINGS THREADS QUERY_THREADS QUERIES REPS COUNTERS
        COUNTER_EVENTS PLACEMENT SSHASH_THREADS SSHASH_M $(compgen -v | grep '^THREADS_'))

mkdir -p "${out_root}"
log_file="${out_root}/feed-${node}.log"
done_file="${out_root}/feed-${node}.done"
touch "${done_file}"


# Logs the message $*.
log()
{
    echo "$(date '+%F %T') $*" | tee -a "${log_file}"
}


# Prints the number of outputs in the directory $1: its runs, distinct counts, and SSHash index costs.
outputs()
{
    find "${1}" -mindepth 2 -maxdepth 2 -name '*.tsv' 2> /dev/null | wc -l
}


# Submits a job of the experiment $1, and sets `id` to its ID. sbatch can fail after queueing the job, so a failed
# submission is retried only if no job of its name shows up within a minute. Returns 1 if submitting fails.
submit()
{
    local experiment="${1}" name="kh-${1}-${node}" since out try

    for try in 1 2 3 4 5; do
        since=$(date +%FT%T)
        out=$(env $(printf -- "-u %s " "${matrix[@]}") EXPERIMENTS="${experiment}" sbatch --parsable "${where[@]}" \
            --nodelist="${node}" --job-name="${name}" --output="${out_root}/${experiment}/bench-%j.out" \
            scripts/bench/explorer.sbatch "${data_dir}" "${out_root}/${experiment}" "${plan}" 2>&1)
        id=$(grep -E -o '^[0-9]+' <<< "${out}" | tail -n 1)
        if [ -z "${id}" ]; then
            log "sbatch of ${name} failed: ${out}"
            sleep 60
            id=$(sacct -n -X -u "${USER}" --name="${name}" -S "${since}" -o JobID | grep -E -o '^[0-9]+' | head -n 1)
        fi

        [ -n "${id}" ] && return 0
    done

    return 1
}


# Waits for the job `id` to end, and sets `state` and `code` to its end state and exit code.
await()
{
    while :; do
        IFS='|' read -r code state < <(sacct -j "${id}" -X -n -P -o ExitCode,State 2> /dev/null | head -n 1)
        case "${state}" in
            ""|PENDING|RUNNING|REQUEUED|SUSPENDED|COMPLETING|CONFIGURING|RESIZING) sleep 60 ;;
            *) return ;;
        esac
    done
}


for experiment in "$@"; do
    if grep -q -x -F "${experiment}" "${done_file}"; then
        log "Skipping ${experiment}: done"
        continue
    fi

    out_dir="${out_root}/${experiment}"
    mkdir -p "${out_dir}"
    busy=0
    while :; do
        before=$(outputs "${out_dir}")
        if ! submit "${experiment}"; then
            log "Stopping: could not submit ${experiment}"
            exit 1
        fi

        log "Submitted ${experiment} as job ${id}"
        await
        new=$(( $(outputs "${out_dir}") - before ))
        log "Job ${id} of ${experiment} ended ${state} (exit ${code}), with ${new} new outputs"
        case "${state}" in
            COMPLETED)
                echo "${experiment}" >> "${done_file}"
                break ;;
            TIMEOUT)
                if (( new == 0 )); then
                    log "Stopping: ${experiment} wrote no output within the time limit"
                    exit 1
                fi ;;
            FAILED)
                if [ "${code}" != "4:0" ] || (( ++busy > busy_tries )); then
                    log "Stopping: ${experiment} failed; see ${out_dir}/bench-${id}.out"
                    exit 1
                fi

                sleep "${busy_wait}" ;;
            *)
                log "Stopping: job ${id} of ${experiment} ended ${state}"
                exit 1 ;;
        esac
    done
done

log "Done"
