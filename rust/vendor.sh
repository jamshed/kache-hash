#!/usr/bin/env bash
# Copies the sources that the crate kache-hash-sys compiles, kache-hash's headers and C interface, into its `vendor/`,
# and kache-hash's license into both crates, for `cargo package` and `cargo publish`. Within kache-hash's tree, the
# crate's build takes the tree's sources instead.
# Usage: rust/vendor.sh
set -eu
here="$(cd "$(dirname "$0")" && pwd)"; root="$(cd "${here}/.." && pwd)"
v="${here}/kache-hash-sys/vendor"
rm -rf "${v}"
mkdir -p "${v}/include" "${v}/capi"
cp -R "${root}/include/kache_hash" "${v}/include/"
cp "${root}/capi/"*.cpp "${root}/capi/"*.hpp "${v}/capi/"
for crate in kache-hash-sys kache-hash; do
    cp "${root}/LICENSE" "${here}/${crate}/LICENSE"
done
echo "vendored $(find "${v}" -type f | wc -l | tr -d ' ') files into ${v}"
