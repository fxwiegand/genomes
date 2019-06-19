extern crate rust_htslib;
extern crate bit_vec;

use rust_htslib::bcf;
use rust_htslib::prelude::*;
use std::{fmt, str};
use std::path::Path;
use rust_htslib::bcf::Read;

#[derive(Serialize, Clone)]
pub struct Variant {
    reference: String,
    alternatives: Vec<String>,
    pos: u32,
    row: i8,
}

pub fn read_indexed_vcf(path: &Path, chrom: u8, from: u32, to: u32) -> Vec<Variant> {
    let mut vcf = rust_htslib::bcf::IndexedReader::from_path(&path).unwrap();

    let chr = chrom.to_string();
    let c = chr.as_bytes();
    let rid = vcf.header().name2rid(c).unwrap();

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

        let var = Variant {
            reference: rfrce,
            alternatives: altve,
            pos: pos,
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
        let c = header.rid2name(rid);

        let mut chr_str = String::new();
        chr_str.push_str(str::from_utf8(c).unwrap());

        let chrom: u8 = chr_str.parse().unwrap();


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

        let var = Variant {
            reference: rfrce,
            alternatives: altve,
            pos: pos,
            row: -1,
        };

        variants.push(var);
    }

    variants
}