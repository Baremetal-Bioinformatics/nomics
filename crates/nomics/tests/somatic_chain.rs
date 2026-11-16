#[test]
fn small_fixture_chain_matches_subcommands() {
    let status = std::process::Command::new("python3")
        .arg(concat!(env!("CARGO_MANIFEST_DIR"), "/../../scripts/somatic_chain_fixture.py"))
        .arg(env!("CARGO_BIN_EXE_nomics"))
        .status().unwrap();
    assert!(status.success());
}
