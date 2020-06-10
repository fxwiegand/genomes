use std::fs;
use std::process::Command;

/// Compare an output file to the expected output and delete the output file.
fn test_output(result: &str, expected: &str) {
    assert!(Command::new("cmp")
        .arg(result)
        .arg(expected)
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success());
    fs::remove_file(result).unwrap();
}

#[test]
fn test_report() {
    assert!(
        Command::new("bash")
            .arg("-c")
            .arg("target/release/genomes report -r tests/resources/test.bam tests/resources/ref.fa tests/resources/report-test.vcf.gz chr1 > tests/report.html")
            .spawn()
            .unwrap()
            .wait()
            .unwrap()
            .success()
    );
    test_output("tests/report.html", "tests/expected/report.html");
}
