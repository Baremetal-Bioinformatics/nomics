//! Exercise policy and append integrity across independent test processes.
#[path = "../../../test-support/fixtures.rs"]
mod fixture;
use std::process::Command;

#[test]
fn fixture_child() {
    if let Ok(mode) = std::env::var("NOMICS_FIXTURE_TEST_CHILD") {
        fixture::check(mode == "executed", "test oracle: provide /example/oracle");
    }
}

#[test]
fn policy_and_concurrent_append() {
    let dir = std::env::temp_dir().join(format!("nomics-fixture-policy-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let log = dir.join("decisions.tsv");
    let _ = std::fs::remove_file(&log);
    let child = |mode: &str, strict: bool| {
        let mut cmd = Command::new(std::env::current_exe().unwrap());
        cmd.args(["--exact", "fixture_child", "--nocapture"])
            .env("NOMICS_FIXTURE_TEST_CHILD", mode)
            .env("NOMICS_FIXTURE_LOG", &log)
            .env_remove("NOMICS_REQUIRE_FIXTURES");
        if strict {
            cmd.env("NOMICS_REQUIRE_FIXTURES", "1");
        }
        cmd
    };
    assert!(child("skipped", false).output().unwrap().status.success());
    let failure = child("skipped", true).output().unwrap();
    assert!(!failure.status.success());
    let message = String::from_utf8_lossy(&failure.stderr);
    assert!(message.contains("required fixture/oracle/binary missing"));
    assert!(message.contains("/example/oracle"));
    let mut processes = Vec::new();
    for _ in 0..8 {
        processes.push(child("executed", true).spawn().unwrap());
    }
    for mut process in processes {
        assert!(process.wait().unwrap().success());
    }
    let records = std::fs::read_to_string(&log).unwrap();
    assert_eq!(records.lines().count(), 10);
    assert_eq!(
        records
            .lines()
            .filter(|s| s.contains("\texecuted\t"))
            .count(),
        8
    );
    assert_eq!(
        records
            .lines()
            .filter(|s| s.contains("\tskipped\t"))
            .count(),
        2
    );
    assert!(records.lines().all(|s| s.split('\t').count() == 4));
    std::fs::remove_dir_all(dir).unwrap();
}
