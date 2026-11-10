#!/usr/bin/env bash
# Checks kache-hash as a dependency. Configures, builds and installs a copy of the source tree, static and then shared,
# each into a prefix of its own, with no network access where the platform can deny it (and the copy must get no
# `external/`, where downloads would go). Then builds and runs the consumers here against each install, through
# `find_package` and through pkg-config, and against the source tree, through FetchContent: from source, and with
# FIND_PACKAGE_ARGS taking the static install.
# Usage: test/consumer/run.sh <work-dir> [CMake arguments for each configure...], e.g. -DCMAKE_OSX_ARCHITECTURES=x86_64
# on Apple silicon, for builds that run under Rosetta.
set -eu
here="$(cd "$(dirname "$0")" && pwd)"; root="$(cd "${here}/../.." && pwd)"
w="$1"; shift
rm -rf "${w}"; mkdir -p "${w}"
src="${w}/src"
rsync -a --exclude '/build*/' --exclude /external/ --exclude /data/ --exclude /local/ --exclude /.git --exclude target/ \
    "${root}/" "${src}/"
j="$(getconf _NPROCESSORS_ONLN)"

# No network for kache-hash's own builds: a new network namespace on Linux, a sandbox profile on macOS.
offline=()
if unshare -rn true 2> /dev/null; then
    offline=(unshare -rn)
elif sandbox-exec -p '(version 1)(allow default)(deny network*)' true 2> /dev/null; then
    offline=(sandbox-exec -p '(version 1)(allow default)(deny network*)')
fi
echo "kache-hash builds offline: ${offline[*]:-not enforced on this platform}"

# Direct compiles for pkg-config take the architecture of the CMake builds.
cflags=""
[[ " $* " == *" -DCMAKE_OSX_ARCHITECTURES=x86_64 "* ]] && cflags="-arch x86_64"

cases=0; fail=0
run() # Runs the consumers `$3...` built in the directory `$2`, as the case `$1`.
{
    local name="$1" dir="$2" b out; shift 2
    for b in "$@"; do
        cases=$((cases + 1))
        if out=$(cd "${dir}" && "./${b}" 2>&1); then
            echo "ok    ${name} ${b}: ${out}"
        else
            echo "FAIL  ${name} ${b}: ${out}"; fail=$((fail + 1))
        fi
    done
}

build() # Configures the consumers into `$1` with the CMake arguments `$2...`, and builds them.
{
    local dir="$1"; shift
    cmake -S "${here}" -B "${dir}" "$@" > "${dir}.log" 2>&1 && cmake --build "${dir}" -j "${j}" >> "${dir}.log" 2>&1 \
        || { echo "FAIL  consumers in ${dir}; see ${dir}.log"; tail -20 "${dir}.log"; fail=$((fail + 1)); return 1; }
}

for type in static shared; do
    shared=OFF; [ "${type}" = shared ] && shared=ON
    p="${w}/prefix-${type}"
    if ! { "${offline[@]}" cmake -S "${src}" -B "${w}/build-${type}" -DBUILD_SHARED_LIBS=${shared} \
                -DCMAKE_INSTALL_PREFIX="${p}" -DKACHE_HASH_BUILD_TESTS=OFF "$@" && "${offline[@]}" cmake --build "${w}/build-${type}" -j "${j}" &&
           "${offline[@]}" cmake --install "${w}/build-${type}"; } > "${w}/build-${type}.log" 2>&1; then
        echo "FAIL  kache-hash ${type} build; see ${w}/build-${type}.log"; tail -20 "${w}/build-${type}.log"; exit 1
    fi
    lib="$(dirname "$(find "${p}" -name 'libkache_capi.*' | head -1)")"
    echo "installed ${type}: $(cd "${p}" && find . -type f -o -type l | sort | grep -v '/include/kache_hash/' | tr '\n' ' ')"

    build "${w}/find-${type}" -DCMAKE_PREFIX_PATH="${p}" "$@" && run "find_package/${type}" "${w}/find-${type}" table capi counter

    export PKG_CONFIG_PATH="${lib}/pkgconfig:${p}/share/pkgconfig"
    mkdir -p "${w}/pkg-config-${type}"
    if cc ${cflags} -std=c99 -Wall -Wextra -Wpedantic -Werror "${here}/capi.c" -o "${w}/pkg-config-${type}/capi" \
            $(pkg-config --cflags --libs kache_capi) -Wl,-rpath,"${lib}" &&
       c++ ${cflags} -std=c++20 -mavx2 "${here}/table.cpp" -o "${w}/pkg-config-${type}/table" \
            $(pkg-config --cflags --libs kache_hash); then
        run "pkg-config/${type}" "${w}/pkg-config-${type}" table capi
    else
        echo "FAIL  pkg-config/${type} builds"; fail=$((fail + 1))
    fi
done

build "${w}/fetch-source" -DKACHE_HASH_SOURCE="${src}" "$@" && run "FetchContent/source" "${w}/fetch-source" table capi counter
build "${w}/fetch-found" -DKACHE_HASH_SOURCE="${src}" -DCMAKE_PREFIX_PATH="${w}/prefix-static" "$@" &&
    run "FetchContent/found" "${w}/fetch-found" table capi counter
found="$(grep '^kache_hash_DIR:' "${w}/fetch-found/CMakeCache.txt" | cut -d= -f2-)"
if [[ "${found}" != "${w}/prefix-static/"* ]]; then
    echo "FAIL  FetchContent/found took kache_hash from ${found:-its source}, not the install"; fail=$((fail + 1))
fi
if [ -e "${src}/external" ]; then
    echo "FAIL  the builds wrote ${src}/external, where downloads go"; fail=$((fail + 1))
fi

echo "consumers: ${cases} runs, failures: ${fail}"
exit $((fail > 0))
