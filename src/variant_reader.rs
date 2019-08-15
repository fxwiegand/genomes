extern crate rust_htslib;
extern crate bit_vec;

use std::{str};
use std::path::Path;
use rust_htslib::bcf::Read;

#[derive(Serialize, Clone)]
pub struct Variant {
    marker_type: String,
    reference: String,
    alternatives: Vec<String>,
    start_position: f32,
    end_position: f32,
    row: i8,
}

pub fn read_indexed_vcf(path: &Path, chrom: String, from: u32, to: u32) -> Vec<Variant> {
    let mut vcf = rust_htslib::bcf::IndexedReader::from_path(&path).unwrap();

    let rid = vcf.header().name2rid(chrom.as_bytes()).unwrap();

    vcf.fetch(rid, from, to).unwrap();

    let mut variants: Vec<Variant> = Vec::new();

    for r in vcf.records() {
        let mut rec = r.unwrap();

        let pos = rec.pos();
        let alleles = rec.alleles();


        let ref_vec = alleles[0].clone();
        let mut rfrce = String::from("");

        for c in ref_vec {
            rfrce.push(*c as char);
        }

        let mut altve: Vec<String> = Vec::new();

        for i in 1..alleles.len() {
            let alt = alleles[i];
            let mut str = String::from("");

            for c in alt {
                str.push(*c as char);
            }

            altve.push(str);
        }

        let var_string = String::from("Variant");

        let var = Variant {
            marker_type: var_string,
            reference: rfrce,
            alternatives: altve,
            start_position: pos as f32 - 0.5,
            end_position: pos as f32 + 0.5,
            row: -1,
        };

        variants.push(var);
    }

    variants
}


pub fn read_vcf(path: &Path) -> Vec<Variant> {
    let mut vcf = rust_htslib::bcf::Reader::from_path(&path).unwrap();

    let mut variants: Vec<Variant> = Vec::new();


    for r in vcf.records() {
        let mut rec = r.unwrap();
        let header = rec.header();
        let rid = rec.rid().unwrap();
        let c = header.rid2name(rid).unwrap();

        let mut chr_str = String::new();
        chr_str.push_str(str::from_utf8(c).unwrap());


        let pos = rec.pos();
        let alleles = rec.alleles();


        let ref_vec = alleles[0].clone();
        let mut rfrce = String::from("");

        for c in ref_vec {
            rfrce.push(*c as char);
        }

        let mut altve: Vec<String> = Vec::new();

        for i in 1..alleles.len() {
            let alt = alleles[i];
            let mut str = String::from("");

            for c in alt {
                str.push(*c as char);
            }

            altve.push(str);
        }

        let var_string = String::from("Variant");

        let var = Variant {
            marker_type: var_string,
            reference: rfrce,
            alternatives: altve,
            start_position: pos as f32 - 0.5,
            end_position: pos as f32 + 0.5,
            row: -1,
        };

        variants.push(var);
    }

    variants
}