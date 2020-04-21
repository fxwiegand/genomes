use super::*;
use std::path::Path;

#[test]
fn reference_test() {
    let ref_bases = read_fasta(Path::new("tests/resources/ref.fa"), String::from("1"), 1, 10);

    let mut compare_ref = Vec::new();

    let bases = String::from("TGCCGGGGT");
    let mut pos: f64 = 1.0;
    for char in bases.chars() {
        let base = Nucleobase {
            start_position: pos - 0.5,
            end_position: pos + 0.5,
            marker_type: char,
            row: 0
        };
        compare_ref.push(base);
        pos += 1.0;
    }
    assert_eq!(compare_ref,ref_bases);
}

#[test]
fn empty_reference_test() {
    let ref_bases = read_fasta(Path::new("tests/resources/ref.fa"), String::from("1"), 11, 11);

    let compare_ref:Vec<Nucleobase> = Vec::new();

    assert_eq!(compare_ref,ref_bases);
}
