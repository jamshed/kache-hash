#!/usr/bin/env bash
# Builds the documentation site into `<out>/site`: the pages of `docs/`, by MkDocs Material, and the API reference
# under `api/`, by Doxygen from the headers. It takes `mkdocs`, from `docs/requirements.txt`, and `doxygen`, on the PATH.
# Usage: docs/build.sh [<out>, default build-docs]
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
out="$(mkdir -p "${1:-${root}/build-docs}" && cd "${1:-${root}/build-docs}" && pwd)"
cd "${root}"

version="$(sed -n -E 's/^#define KACHE_HASH_VERSION_(MAJOR|MINOR|PATCH) ([0-9]+)$/\2/p' include/kache_hash/version.hpp | paste -sd. -)"
rm -rf "${out}/doxygen" "${out}/site"
KACHE_HASH_VERSION="${version}" KACHE_HASH_DOXYGEN_OUTPUT="${out}/doxygen" doxygen docs/api/Doxyfile
mkdocs build --quiet --site-dir "${out}/site"
cp -R "${out}/doxygen/html" "${out}/site/api"
echo "built the site of kache-hash ${version} into ${out}/site"
