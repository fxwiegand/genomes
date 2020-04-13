use super::*;
use std::path::Path;

#[test]
fn duplicate_test() {
    let mut variants = read_indexed_vcf(Path::new("tests/resources/duplicate.vcf.gz"), String::from("11"), 150000, 151000);
    let var = variants.pop().unwrap();

    let mut alt = Vec::new();
    let allel = String::from("ATCATC");
    alt.push(allel);

    let test_variant = Variant {
        marker_type: String::from("Variant"),
        reference: String::from("ATC"),
        alternatives: Some(alt),
        start_position: 150188.5 as f64, // - 1 due to 0-basing, - 0.5 due to start pos
        end_position: 150191.5 as f64, // -1 due to 0-basing, + 0.5 du to end pos
    };
    assert_eq!(var, test_variant);
}


#[test]
fn inversion_test() {
    let mut variants = read_indexed_vcf(Path::new("tests/resources/inversion.vcf.gz"), String::from("11"), 150000, 151000);
    let var = variants.pop().unwrap();

    let mut alt = Vec::new();
    let allel = String::from("CTTAG");
    alt.push(allel);

    let test_variant = Variant {
        marker_type: String::from("Variant"),
        reference: String::from("GATTC"),
        alternatives: Some(alt),
        start_position: 150186.5 as f64, // - 1 due to 0-basing, - 0.5 due to start pos
        end_position: 150191.5 as f64, // -1 due to 0-basing, + 0.5 du to end pos
    };
    assert_eq!(var, test_variant);
}


#[test]
fn deletion_test() {
    let mut variants = read_indexed_vcf(Path::new("tests/resources/deletion.vcf.gz"), String::from("11"), 150000, 151000);
    let var = variants.pop().unwrap();

    let test_variant = Variant {
        marker_type: String::from("Variant"),
        reference: String::from("G"),
        alternatives: None,
        start_position: 150186.5 as f64, // - 1 due to 0-basing, - 0.5 due to start pos
        end_position: 150773.5 as f64, // -1 due to 0-basing, + 0.5 du to end pos
    };
    assert_eq!(var, test_variant);
}