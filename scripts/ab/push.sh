#!/usr/bin/env bash
# Pushes a source tree to explorer's ~/kh-perf/src/<name>/, for A/B jobs (mkjob.sh) and counter jobs (cnt.sh) to build:
# the tree of a commit, or the working tree (WT, the default).
# Usage: scripts/ab/push.sh <name> [<commit>|WT]
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
host=explorer-02
name="$1"; rev="${2:-WT}"
tmp="$(mktemp -d)"
if [ "${rev}" = WT ]; then
    rsync -a --exclude "build*/" --exclude external/ --exclude data --exclude .git --exclude local/ "${root}/" "${tmp}/"
else
    git -C "${root}" archive "${rev}" | tar -x -C "${tmp}"
fi
COPYFILE_DISABLE=1 tar --no-xattrs -C "${tmp}" -czf - . | ssh -o BatchMode=yes ${host} "rm -rf ~/kh-perf/src/${name} && mkdir -p ~/kh-perf/src/${name} && tar --warning=no-unknown-keyword -xzf - -C ~/kh-perf/src/${name}"
rm -rf "${tmp}"
echo "pushed ${name} (${rev})"
