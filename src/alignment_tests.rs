use super::*;
use alignment_reader::Marker;
use std::path::Path;

#[test]
fn match_test() {
    let (_bam, mut matches) = get_reads(
        Path::new("tests/resources/test.bam"),
        Path::new("tests/resources/ref.fa"),
        String::from("chr1"),
        0,
        100,
    );
    matches.retain(|m| m.marker_type == Marker::Match);

    let mut compare_matches = Vec::new();

    let m1 = AlignmentMatch {
        marker_type: Marker::Match,
        start_position: 3.5,
        end_position: 19.5,
        flags: vec![1, 2, 32, 64],
        name: "sim_Som1-5-2_chr1_1_1acd6f".to_string(),
        read_start: 4,
        read_end: 789364,
    };

    compare_matches.push(m1);

    let m2 = AlignmentMatch {
        marker_type: Marker::Match,
        start_position: 19.5,
        end_position: 99.5,
        flags: vec![1, 2, 32, 64],
        name: "sim_Som1-5-2_chr1_1_1acd6f".to_string(),
        read_start: 4,
        read_end: 789364,
    };

    compare_matches.push(m2);

    assert_eq!(compare_matches, matches);
}

#[test]
fn mismatch_test() {
    let (mut bam, _matches) = get_reads(
        Path::new("tests/resources/test.bam"),
        Path::new("tests/resources/ref.fa"),
        String::from("chr1"),
        0,
        110,
    );
    bam.retain(|m| m.marker_type == Marker::T);

    let mut compare_bam = Vec::new();

    let m = AlignmentNucleobase {
        marker_type: Marker::T,
        bases: "T".to_string(),
        start_position: 99.5,
        end_position: 100.5,
        flags: vec![1, 2, 32, 64],
        name: "sim_Som1-5-2_chr1_1_1acd6f".to_string(),
        read_start: 4,
        read_end: 789364,
    };

    compare_bam.push(m);

    assert_eq!(compare_bam, bam);
}

#[test]
fn insertion_test() {
    let (mut bam, _matches) = get_reads(
        Path::new("tests/resources/test.bam"),
        Path::new("tests/resources/ref.fa"),
        String::from("chr1"),
        0,
        100,
    );
    bam.retain(|m| m.marker_type == Marker::Insertion);

    let mut compare_bam = Vec::new();

    let m = AlignmentNucleobase {
        marker_type: Marker::Insertion,
        bases: "AA".to_string(),
        start_position: 19.0,
        end_position: 20.0,
        flags: vec![1, 2, 32, 64],
        name: "sim_Som1-5-2_chr1_1_1acd6f".to_string(),
        read_start: 4,
        read_end: 789364,
    };

    compare_bam.push(m);

    assert_eq!(compare_bam, bam);
}

#[test]
fn deletion_test() {
    let (mut bam, _matches) = get_reads(
        Path::new("tests/resources/del.bam"),
        Path::new("tests/resources/ref.fa"),
        String::from("chr1"),
        0,
        100,
    );
    bam.retain(|m| m.marker_type == Marker::Deletion);

    let mut compare_bam = Vec::new();

    let m = AlignmentNucleobase {
        marker_type: Marker::Deletion,
        bases: "".to_string(),
        start_position: 18.5,
        end_position: 19.5,
        flags: vec![2, 32, 64],
        name: "sim_Som1-5-2_chr1_1_1acd6f".to_string(),
        read_start: 4,
        read_end: 103,
    };

    compare_bam.push(m);

    assert_eq!(compare_bam, bam);
}
