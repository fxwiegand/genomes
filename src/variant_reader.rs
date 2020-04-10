extern crate rust_htslib;
extern crate bit_vec;

use std::path::Path;
use rust_htslib::bcf::Read;
use regex::Regex;

#[derive(Serialize, Clone)]
pub struct Variant {
    pub(crate) marker_type: String,
    pub(crate) reference: String,
    pub(crate) alternatives: Option<Vec<String>>,
    pub(crate) start_position: f64,
    pub(crate) end_position: f64,
}

pub fn read_indexed_vcf(path: &Path, chrom: String, from: u64, to: u64) -> Vec<Variant> {
    let mut vcf = rust_htslib::bcf::IndexedReader::from_path(&path).unwrap();

    let rid = vcf.header().name2rid(chrom.as_bytes()).unwrap();

    vcf.fetch(rid, from, to).unwrap();

    let mut variants: Vec<Variant> = Vec::new();

    for r in vcf.records() {
        let mut rec = r.unwrap();

        let pos = rec.pos();
        let end_pos = match rec.info(b"END").integer() {
            Ok(Some(end_pos)) => {
                let end_pos = end_pos[0] as f64 + 0.5;
                Some(end_pos)
            }
            _ => None,
        };
        let alleles = rec.alleles();


        let ref_vec = alleles[0].clone();
        let mut rfrce = String::from("");

        let mut len: u8 = 0;

        for c in ref_vec {
            rfrce.push(*c as char);
            len += 1;
        }

        if alleles[1] == b"<DEL>" {
            let var_string = String::from("Variant");

            let var = Variant {
                marker_type: var_string,
                reference: rfrce.clone(),
                alternatives: None,
                start_position: pos as f64 - 0.5,
                end_position: end_pos.unwrap(),
            };

            variants.push(var);
        } else {
            let mut altve: Vec<String> = Vec::new();

            for i in 1..alleles.len() {
                let alt = alleles[i];

                let mut str = String::from("");

                for c in alt {
                    str.push(*c as char);
                }

                let cnv = Regex::new(r"^<CN\d>$").unwrap();

                if cnv.is_match(str.as_ref()) {
                    warn!("Use of unsupported Copy-Number-Variation {}", str)
                } else {
                    altve.push(str);
                }
            }

            if altve.len() > 0 {
                let var_string = String::from("Variant");

                let var = Variant {
                    marker_type: var_string,
                    reference: rfrce,
                    alternatives: Some(altve),
                    start_position: pos as f64 - 0.5,
                    end_position: pos as f64 - 0.5 + len as f64,
                };

                variants.push(var);
            }
        }

    }

    variants
}