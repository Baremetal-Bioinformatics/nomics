//! Compiles `shim/*.cu` into a static library when the `cuda` feature is on,
//! or `shim/metal/` (Objective-C shim + embedded metallib) with `metal`.
//! Env overrides: NOMICS_NVCC (default /usr/local/cuda/bin/nvcc),
//! NOMICS_CUDA_ARCH (default sm_121, GB10), NOMICS_METALLIB (a prebuilt
//! metallib of the same source, for hosts without the Metal compiler).
use std::{
    env,
    path::{Path, PathBuf},
    process::Command,
};

const SOURCES: &[&str] = &[
    "shim/smoke.cu",
    "shim/stream.cu",
    "shim/seed.cu",
    "shim/pairhmm.cu",
    "shim/ksw.cu",
    "shim/extend.cu",
    "shim/fused.cu",
];

fn main() {
    println!("cargo:rerun-if-changed=shim");
    println!("cargo:rerun-if-env-changed=NOMICS_NVCC");
    println!("cargo:rerun-if-env-changed=NOMICS_CUDA_ARCH");
    println!("cargo:rerun-if-env-changed=NOMICS_METALLIB");
    let cuda = env::var_os("CARGO_FEATURE_CUDA").is_some();
    let metal = env::var_os("CARGO_FEATURE_METAL").is_some();
    assert!(
        !(cuda && metal),
        "features `cuda` and `metal` are exclusive"
    );
    if metal {
        build_metal();
    }
    if !cuda {
        return;
    }
    let nvcc = env::var("NOMICS_NVCC").unwrap_or_else(|_| "/usr/local/cuda/bin/nvcc".into());
    let arch = env::var("NOMICS_CUDA_ARCH").unwrap_or_else(|_| "sm_121".into());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let lib = out.join("libnomics_gpu_shim.a");
    let status = Command::new(&nvcc)
        .args([
            "-std=c++17",
            "-O3",
            "-lineinfo",
            "--fmad=false",
            "-lib",
            "-Xcompiler",
            "-fPIC",
        ])
        .arg(format!("-arch={arch}"))
        .arg("-I/usr/local/cuda/include/cccl")
        .arg("-o")
        .arg(&lib)
        .args(SOURCES)
        .status()
        .unwrap_or_else(|e| panic!("failed to run {nvcc}: {e}; set NOMICS_NVCC"));
    assert!(status.success(), "nvcc failed");
    let cuda_root = PathBuf::from(&nvcc)
        .parent()
        .and_then(|p| p.parent())
        .unwrap()
        .to_path_buf();
    println!("cargo:rustc-link-search=native={}", out.display());
    println!(
        "cargo:rustc-link-search=native={}",
        cuda_root.join("lib64").display()
    );
    println!("cargo:rustc-link-lib=static=nomics_gpu_shim");
    println!("cargo:rustc-link-lib=dylib=cudart");
    println!("cargo:rustc-link-lib=dylib=stdc++");
}

const METAL_KERNELS: &[&str] = &[
    "shim/metal/smoke.metal",
    "shim/metal/ksw.metal",
    "shim/metal/seed.metal",
    "shim/metal/fused.metal",
];

/// `xcrun <args>`, panicking with the command line on failure.
fn xcrun(args: &[&str]) {
    let status = Command::new("xcrun")
        .args(args)
        .status()
        .unwrap_or_else(|e| panic!("failed to run xcrun: {e}"));
    assert!(status.success(), "xcrun {} failed", args.join(" "));
}

/// METAL_KERNELS → one metallib at `metallib`.
fn compile_kernels(out: &Path, metallib: &Path) {
    let mut airs = Vec::new();
    for src in METAL_KERNELS {
        let air = out.join(Path::new(src).with_extension("air").file_name().unwrap());
        // Fast math is off: kernels must match the CPU port bit for bit.
        xcrun(&[
            "-sdk",
            "macosx",
            "metal",
            "-std=metal3.1",
            // The library records its minimum OS; without this it is the
            // build host's, and an older macOS refuses a prebuilt metallib.
            "-mmacosx-version-min=15.0",
            "-O3",
            "-fno-fast-math",
            "-c",
            src,
            "-o",
            air.to_str().unwrap(),
        ]);
        airs.push(air.to_str().unwrap().to_owned());
    }
    let mut args = vec!["-sdk", "macosx", "metallib", "-o", metallib.to_str().unwrap()];
    args.extend(airs.iter().map(String::as_str));
    xcrun(&args);
}

/// The kernels' metallib in OUT_DIR: compiled from METAL_KERNELS, or copied
/// from NOMICS_METALLIB (built from the same source elsewhere: the Metal
/// compiler is a separate download that a host may lack).
fn metallib(out: &Path) -> PathBuf {
    let metallib = out.join("nomics.metallib");
    if let Some(prebuilt) = env::var_os("NOMICS_METALLIB") {
        println!("cargo:rerun-if-changed={}", Path::new(&prebuilt).display());
        std::fs::copy(&prebuilt, &metallib)
            .unwrap_or_else(|e| panic!("NOMICS_METALLIB {prebuilt:?}: {e}"));
    } else {
        compile_kernels(out, &metallib);
    }
    metallib
}

/// The metallib, embedded as a C byte array next to the ObjC shim, so the
/// binary carries its kernels (no runtime path lookup).
fn build_metal() {
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let path = |p: &PathBuf| p.to_str().unwrap().to_owned();
    let metallib = metallib(&out);
    let bytes = std::fs::read(&metallib).unwrap();
    let mut c = String::from("#include <stddef.h>\nconst unsigned char nomics_metallib[] = {");
    for (i, b) in bytes.iter().enumerate() {
        if i % 24 == 0 {
            c.push('\n');
        }
        c.push_str(&format!("{b},"));
    }
    c.push_str(&format!(
        "\n}};\nconst size_t nomics_metallib_len = {};\n",
        bytes.len()
    ));
    let data_c = out.join("metallib_data.c");
    std::fs::write(&data_c, c).unwrap();

    let shim_o = out.join("nomics_metal.o");
    let data_o = out.join("metallib_data.o");
    xcrun(&[
        "clang",
        "-c",
        "-O2",
        "-fobjc-arc",
        "-Wall",
        "-Werror",
        "shim/metal/nomics_metal.m",
        "-o",
        &path(&shim_o),
    ]);
    xcrun(&["clang", "-c", "-O2", &path(&data_c), "-o", &path(&data_o)]);
    let lib = out.join("libnomics_gpu_shim.a");
    let _ = std::fs::remove_file(&lib);
    xcrun(&[
        "libtool",
        "-static",
        "-o",
        &path(&lib),
        &path(&shim_o),
        &path(&data_o),
    ]);

    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=nomics_gpu_shim");
    for framework in ["Metal", "Foundation", "IOKit", "CoreFoundation"] {
        println!("cargo:rustc-link-lib=framework={framework}");
    }
}
