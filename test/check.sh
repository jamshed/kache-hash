#!/usr/bin/env bash
# Builds kache-hash's tests with CMake in build-local/check, as x86-64 (under Rosetta on Apple silicon), and runs the
# quick ones, on the small inputs in test/data/; with `full`, the ones on the datasets of `scripts/fetch-data.sh small`
# in data/ too; with `consumer`, the consumers of kache-hash's installs and sources instead. Extra arguments go to
# CMake's configuration, e.g. -DKACHE_HASH_TEST_THREADS=4.
# Usage: test/check.sh [quick|full|consumer] [CMake arguments...]
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
labels='^quick$'
case "${1:-}" in
    quick) shift ;;
    full) labels='^(quick|full)$'; shift ;;
    consumer) labels='^consumer$'; shift ;;
esac
b="${root}/build-local/check"
args=()
[ "$(uname -s)" = Darwin ] && [ "$(uname -m)" = arm64 ] && args+=(-DCMAKE_OSX_ARCHITECTURES=x86_64)
cpus="$(getconf _NPROCESSORS_ONLN)"
mkdir -p "${b}"
cmake --log-level=WARNING -S "${root}" -B "${b}" -DKACHE_HASH_BUILD_TESTS=ON ${args[@]+"${args[@]}"} "$@"
cmake --build "${b}" -j "${cpus}" > "${b}/build.log" 2>&1 || { tail -40 "${b}/build.log"; exit 1; }
ctest --test-dir "${b}" -L "${labels}" -j "${cpus}" --output-on-failure
