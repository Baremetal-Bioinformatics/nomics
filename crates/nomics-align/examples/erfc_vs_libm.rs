//! Compare nomics_align::libm::erfc with libm samples from
//! reference/oracle/erfc_samples.c. usage: erfc_vs_libm <samples.bin>
fn main() {
    let path = std::env::args().nth(1).expect("samples.bin");
    let data = std::fs::read(path).unwrap();
    let (mut n, mut bad) = (0u64, 0u64);
    for rec in data.chunks_exact(16) {
        let x = f64::from_le_bytes(rec[..8].try_into().unwrap());
        let y = f64::from_le_bytes(rec[8..].try_into().unwrap());
        let ours = nomics_align::libm::erfc(x);
        n += 1;
        if ours.to_bits() != y.to_bits() {
            if bad < 10 {
                eprintln!("x={x:e} ({:#x}): libm {y:e} ours {ours:e}", x.to_bits());
            }
            bad += 1;
        }
    }
    println!("erfc: {n} samples, {bad} bit mismatches");
    std::process::exit((bad != 0) as i32);
}
