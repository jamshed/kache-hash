// Compiles kache-hash's C interface library, `kache_capi`, with the `cc` crate, and links it with the C++ standard
// library: from kache-hash's tree when the crate builds within it, and otherwise from the copy of the sources that
// `rust/vendor.sh` puts in the package. The crate's docs say which ISA it targets.

use std::env;
use std::path::{Path, PathBuf};

/// The ISAs that `KACHE_HASH_ARCH` names, as kache-hash's CMake does.
const ARCHES: [&str; 4] = ["native", "x86-64-v3", "x86-64-v4", "none"];

/// The AVX-512 features of the x86-64-v4 level.
const AVX512: [&str; 5] = ["avx512f", "avx512bw", "avx512cd", "avx512dq", "avx512vl"];

fn main() {
    println!("cargo:rerun-if-env-changed=KACHE_HASH_ARCH");

    // docs.rs documents the crate, which takes no linking.
    if env::var_os("DOCS_RS").is_some() {
        return;
    }

    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap();
    if target_arch != "x86_64" {
        panic!("kache-hash needs an x86-64 CPU with AVX2 for now, and this build targets {target_arch}.");
    }

    let root = sources();
    let arch = arch();

    let mut build = cc::Build::new();
    build
        .cpp(true)
        .std("c++20")
        .include(root.join("include"))
        // Release's flags, whatever Cargo's profile: others would compile the table without optimizations, or with its
        // assertions.
        .opt_level(3)
        .debug(false)
        .define("NDEBUG", None)
        .flag("-funroll-loops")
        // Only the interface's symbols are exported, as from kache-hash's CMake.
        .flag("-fvisibility=hidden")
        .flag("-fvisibility-inlines-hidden");

    // An x86-64 build on another host, such as one on Apple silicon to run under Rosetta, can't target the native CPU;
    // it targets Haswell's AVX2 then, as kache-hash's CMake does.
    match arch.as_str() {
        "native" if !env::var("HOST").unwrap().starts_with("x86_64") => build.flag("-march=haswell"),
        "native" => build.flag("-march=native"),
        "none" => &mut build,
        level => build.flag(format!("-march={level}")),
    };

    // Keep branches from crossing or ending on 32-byte boundaries, for Skylake-derived cores' JCC erratum, as kache-hash's
    // CMake does.
    if arch != "none" {
        let pad = if build.get_compiler().is_like_clang() {
            "-mbranches-within-32B-boundaries"
        } else {
            "-Wa,-mbranches-within-32B-boundaries"
        };
        build.flag_if_supported(pad);
    }

    let capi = root.join("capi");
    let mut files: Vec<PathBuf> = capi
        .read_dir()
        .expect("kache-hash's C interface sources")
        .map(|e| e.expect("a C interface source").path())
        .filter(|p| p.extension().is_some_and(|x| x == "cpp"))
        .collect();
    files.sort();
    build.files(&files).compile("kache_capi");

    // The header `kache_hash/kache_capi.h` is under `DEP_KACHE_CAPI_INCLUDE`, for dependents' build scripts.
    println!("cargo:include={}", root.join("include").display());
    for p in [root.join("include").join("kache_hash"), capi] {
        println!("cargo:rerun-if-changed={}", p.display());
    }
}

/// Returns the root of the sources to compile: kache-hash's tree, if the crate is its `rust/kache-hash-sys`, and the
/// package's `vendor` directory otherwise.
fn sources() -> PathBuf {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let tree = manifest.join("../..");
    let within = |root: &Path| {
        root.join("include/kache_hash/kache_capi.h").is_file()
            && root.join("rust/kache-hash-sys").canonicalize().ok() == manifest.canonicalize().ok()
    };
    if within(&tree) {
        return tree;
    }

    let vendor = manifest.join("vendor");
    if !vendor.join("include/kache_hash/kache_capi.h").is_file() {
        panic!("kache-hash-sys has no sources: the package lacks `vendor/`, which `rust/vendor.sh` makes.");
    }
    vendor
}

/// Returns the ISA to compile for: `KACHE_HASH_ARCH`'s, if set; otherwise the native CPU if rustc targets it, x86-64-v4 if
/// rustc's target features include AVX-512's, and x86-64-v3 otherwise, the least that kache-hash supports.
fn arch() -> String {
    if let Ok(a) = env::var("KACHE_HASH_ARCH") {
        if !ARCHES.contains(&a.as_str()) {
            panic!("KACHE_HASH_ARCH is {a}, not one of {}.", ARCHES.join(", "));
        }
        return a;
    }

    if rustc_target_cpu().as_deref() == Some("native") {
        return "native".to_owned();
    }

    let features = env::var("CARGO_CFG_TARGET_FEATURE").unwrap_or_default();
    let features: Vec<&str> = features.split(',').collect();
    let level = if AVX512.iter().all(|f| features.contains(f)) { "x86-64-v4" } else { "x86-64-v3" };
    level.to_owned()
}

/// Returns the CPU that rustc targets, by the last `-C target-cpu` of its flags, if any.
fn rustc_target_cpu() -> Option<String> {
    let flags = env::var("CARGO_ENCODED_RUSTFLAGS").unwrap_or_default();
    let mut args = flags.split('\x1f');
    let mut cpu = None;
    while let Some(a) = args.next() {
        let option = match a {
            "-C" | "--codegen" => args.next(),
            _ => a.strip_prefix("-C").or_else(|| a.strip_prefix("--codegen=")),
        };
        if let Some(c) = option.and_then(|o| o.strip_prefix("target-cpu=")) {
            cpu = Some(c.to_owned());
        }
    }
    cpu
}
