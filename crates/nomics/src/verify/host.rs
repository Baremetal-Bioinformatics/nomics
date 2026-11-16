//! What the report says about the machine: CPU, GPU, memory, OS and kernel.

use std::process::Command;

#[derive(Debug, Clone)]
pub struct Host {
    pub arch: String,
    pub cpu: String,
    pub cores: usize,
    pub memory_gib: u64,
    pub gpu: Option<String>,
    pub os: String,
    pub kernel: String,
}

fn read(path: &str) -> String {
    std::fs::read_to_string(path).unwrap_or_default()
}

fn cpu_model() -> String {
    let info = read("/proc/cpuinfo");
    for key in ["model name", "Hardware", "Model"] {
        if let Some(v) = info
            .lines()
            .find(|l| l.starts_with(key))
            .and_then(|l| l.split_once(':'))
            .map(|(_, v)| v.trim().to_owned())
            .filter(|v| !v.is_empty())
        {
            return v;
        }
    }
    // arm64 kernels often name no model; lscpu decodes the part number.
    Command::new("lscpu")
        .output()
        .ok()
        .and_then(|o| {
            String::from_utf8_lossy(&o.stdout)
                .lines()
                .find(|l| l.starts_with("Model name:"))
                .map(|l| l["Model name:".len()..].trim().to_owned())
        })
        .unwrap_or_else(|| "unknown".to_owned())
}

fn memory_gib() -> u64 {
    read("/proc/meminfo")
        .lines()
        .find_map(|l| l.strip_prefix("MemTotal:"))
        .and_then(|v| v.split_whitespace().next()?.parse::<u64>().ok())
        .map_or(0, |kib| kib / (1 << 20))
}

fn os_name() -> String {
    read("/etc/os-release")
        .lines()
        .find_map(|l| l.strip_prefix("PRETTY_NAME="))
        .map(|v| v.trim_matches('"').to_owned())
        .unwrap_or_else(|| std::env::consts::OS.to_owned())
}

fn gpu_name() -> Option<String> {
    nomics_gpu::device_info(0).ok().map(|i| i.name)
}

pub fn probe() -> Host {
    Host {
        arch: std::env::consts::ARCH.to_owned(),
        cpu: cpu_model(),
        cores: std::thread::available_parallelism().map_or(1, |n| n.get()),
        memory_gib: memory_gib(),
        gpu: gpu_name(),
        os: os_name(),
        kernel: read("/proc/sys/kernel/osrelease").trim().to_owned(),
    }
}
