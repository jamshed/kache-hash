#!/bin/bash

# Runs the measurement driver once per event group of an events file, keeping each run's output as
# `<out-dir>/<group>.tsv`.

events="${1}"
out_dir="${2}"

if [ $# -lt 3 ]; then
    echo "Usage: $0 <events-file> <out-dir> <driver> <driver-args>..." >&2
    exit 1
fi

shift 2
mkdir -p "${out_dir}"

while read -r group spec; do
    echo "Running group ${group}"
    if ! "$@" --perf "${spec}" > "${out_dir}/${group}.tsv" < /dev/null; then
        echo "Failed running group ${group}" >&2
        exit 1
    fi
done < <(grep -vE '^[[:space:]]*(#|$)' "${events}")
