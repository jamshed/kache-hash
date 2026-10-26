#!/usr/bin/env bash
# Creates and submits an A/B job on explorer, which runs ab-run.sh on one node.
# Usage: scripts/ab/mkjob.sh <job-name> <variants-file> <configs-file> <params> [sbatch option...]
# The files and the params are as ab-run.sh describes. Each extra option, e.g. `--nodelist=d0007`, becomes an `#SBATCH`
# line. The environment overrides the defaults for Cascade Lake nodes in `sharing`: KH_PARTITION, KH_TIME, KH_CPUS,
# KH_CONSTRAINT and KH_EXCLUDE (empty to omit). KH_PRE is a command to run in the job first.
# The job directory is ~/kh-perf/jobs/<job-name>/; its Slurm log is job-<id>.out there. Submission goes through
# ~/kache-hash-jobs/submit.sh, which retries sbatch's spurious failures.
set -eu
host=explorer-02
name="$1"; variants="$2"; configs="$3"; params="$4"; shift 4
jd="kh-perf/jobs/${name}"
ssh -o BatchMode=yes ${host} "mkdir -p ~/${jd}"
scp -q -o BatchMode=yes "${variants}" "${host}:${jd}/variants"
scp -q -o BatchMode=yes "${configs}" "${host}:${jd}/configs"
echo "${params}" | ssh -o BatchMode=yes ${host} "cat > ~/${jd}/params"
scp -q -o BatchMode=yes "$(dirname "$0")/ab-run.sh" "${host}:${jd}/"   # A per-job copy: bash reads a running script as it goes.
partition="${KH_PARTITION-sharing}"; time="${KH_TIME-00:59:00}"; cpus="${KH_CPUS-56}"
constraint="${KH_CONSTRAINT-cascadelake}"; exclude="${KH_EXCLUDE-d3230,d3231,d4067}"
extra=""
[ -n "${constraint}" ] && extra+="#SBATCH --constraint=${constraint}"$'\n'
[ -n "${exclude}" ] && extra+="#SBATCH --exclude=${exclude}"$'\n'
for o in "$@"; do extra+="#SBATCH ${o}"$'\n'; done
pre=""
[ -n "${KH_PRE-}" ] && pre="${KH_PRE}"$'\n'
cat <<EOS | ssh -o BatchMode=yes ${host} "cat > ~/${jd}/job.sbatch"
#!/bin/bash
#SBATCH --job-name=khp-${name}
#SBATCH --partition=${partition}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=${time}
#SBATCH --output=job-%j.out
${extra}${pre}bash ~/${jd}/ab-run.sh ~/${jd}
EOS
ssh -o BatchMode=yes ${host} "cd ~/${jd} && ~/kache-hash-jobs/submit.sh job.sbatch"
