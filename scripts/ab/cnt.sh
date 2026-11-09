#!/usr/bin/env bash
# Counts the per-phase events of several pushed variants (push.sh) of the measure driver in a small Slurm job on one
# Cascade Lake core, and records where their instructions retire, per source line and per symbol. Instruction counts per
# phase are deterministic, so they compare across nodes; cycles vary ~10% on shared nodes. cnt-sum.py summarizes.
# Usage: scripts/ab/cnt.sh <tag> "<driver args>" <name>:<source>[:<flags>]...
# Flags within a variant are separated by commas, e.g. a:base:-DFOO,-DBAR=1; '@' stands for a comma within a flag.
# KH_CNT_GROUPS overrides the event groups (space-separated); KH_CNT_REC=0 skips the instruction profiles, and
# KH_CNT_REC_EVENT and KH_CNT_REC_PERIOD set the profiles' sampled event (instructions:u) and period (200003);
# KH_CNT_REC_OPTS adds options, e.g. `-D 18000` to start recording 18 s in.
# The job directory is ~/kh-perf/cnt/<tag>/.
set -eu
host=explorer-02
tag="$1"; args="$2"; shift 2
jd="kh-perf/cnt/${tag}"
groups="${KH_CNT_GROUPS:-br_inst_retired.all_branches=r00c4,br_misp_retired.all_branches=r00c5,uops_issued.any=r010e mem_inst_retired.all_loads=r81d0,mem_inst_retired.all_stores=r82d0,mem_load_retired.l3_miss=r20d1}"
rec="${KH_CNT_REC:-1}"
ssh -o BatchMode=yes ${host} "rm -rf ~/${jd} && mkdir -p ~/${jd}"
{
cat <<EOS
#!/bin/bash
#SBATCH --job-name=khct-${tag}
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --constraint=cascadelake
#SBATCH --exclude=d3230,d3231,d4067
#SBATCH --output=job.out
cd ~/${jd}
lscpu | grep -q "Platinum 8276" || { echo "Not a Cascade Lake node: \$(hostname)"; exit 3; }
hostname
cpu=\$(taskset -pc \$\$ | awk -F': ' '{ print \$2 }' | cut -d, -f1 | cut -d- -f1)
EOS
for v in "$@"; do
    name="${v%%:*}"; rest="${v#*:}"; src="${rest%%:*}"; flags=""; [ "${rest}" != "${src}" ] && flags="$(echo "${rest#*:}" | tr ',' ' ' | tr '@' ',')"
    # The definitions are older trees': the headers set their own since v0.8.0.
    cat <<EOS
rsync -a --exclude "build*/" --exclude external/ ~/kh-perf/src/${src}/ src-${name}/
g++ -std=c++20 -O3 -DNDEBUG -funroll-loops -march=native -Wa,-mbranches-within-32B-boundaries -g -DXXH_INLINE_ALL -DXXH_STATIC_LINKING_ONLY -DXXH_FORCE_MEMORY_ACCESS=0 -DDEVELOP -DL1_CACHE_LINE_SIZE=64 ${flags} -I src-${name}/include -I src-${name}/support -I src-${name}/benchmark src-${name}/benchmark/measure.cpp -o ${name}.bin ~/kh-perf/jemalloc/libjemalloc.a -ldl -lz -pthread || echo "BUILD FAILED: ${name}"
EOS
done
for v in "$@"; do
    name="${v%%:*}"
    i=0
    for g in ${groups}; do
        echo "taskset -c \${cpu} ./${name}.bin ${args} --perf ${g} > ${name}.g${i}.tsv 2> ${name}.g${i}.err; echo \"${name} g${i}: \$(grep -E 'mops' ${name}.g${i}.tsv | tr '\n' ' ')\""
        i=$((i + 1))
    done
    [ "${rec}" = 1 ] || continue
    echo "perf record ${KH_CNT_REC_OPTS:-} -e ${KH_CNT_REC_EVENT:-instructions:u} -c ${KH_CNT_REC_PERIOD:-200003} -o ${name}.perf.data taskset -c \${cpu} ./${name}.bin ${args} > ${name}.rec.tsv 2> ${name}.rec.err"
    echo "perf report -i ${name}.perf.data --stdio --no-children --sort srcline 2>/dev/null | grep -v '^#' | grep -v '^\$' | head -80 > ${name}.srcline.txt"
    echo "perf report -i ${name}.perf.data --stdio --no-children --sort sym 2>/dev/null | grep -v '^#' | grep -v '^\$' | head -30 > ${name}.sym.txt"
done
echo "echo done"
} | ssh -o BatchMode=yes ${host} "cat > ~/${jd}/job.sbatch"
ssh -o BatchMode=yes ${host} "cd ~/${jd} && sbatch job.sbatch"
