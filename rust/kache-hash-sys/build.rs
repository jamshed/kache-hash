// Builds kache-hash's C interface library with kache-hash's own CMake, and links it, with the C++ standard library.

use std::env;
use std::path::PathBuf;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let root = manifest.join("../..").canonicalize().expect("kache-hash's source tree");

    // Release, whatever Cargo's profile: other build types would compile the table without optimizations, or with
    // its assertions.
    let mut cfg = cmake::Config::new(manifest.join("cmake"));
    cfg.define("KACHE_HASH_ROOT", &root).profile("Release");

    // An x86-64 build on Apple silicon runs under Rosetta; kache-hash's CMake then targets Haswell's AVX2.
    let target = env::var("TARGET").unwrap();
    if target == "x86_64-apple-darwin" {
        cfg.define("CMAKE_OSX_ARCHITECTURES", "x86_64");
    }

    let dst = cfg.build();
    println!("cargo:rustc-link-search=native={}", dst.join("lib").display());
    println!("cargo:rustc-link-lib=static=kache_capi");
    if target.contains("apple") {
        println!("cargo:rustc-link-lib=dylib=c++");
    } else {
        println!("cargo:rustc-link-lib=dylib=stdc++");
    }

    for p in ["include", "capi", "CMakeLists.txt"] {
        println!("cargo:rerun-if-changed={}", root.join(p).display());
    }
    println!("cargo:rerun-if-changed={}", manifest.join("cmake").display());
}
