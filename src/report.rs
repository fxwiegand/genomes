use std::path::Path;
use std::error::Error;
use rust_htslib::bcf::Read;
use rustc_serialize::json::Json;
use static_reader::{StaticVariant, get_static_reads};
use variant_reader::{VariantType};
use json_generator::manipulate_json;
use fasta_reader::read_fasta;


#[derive(Serialize, Clone, Debug, PartialEq)]
pub struct Report {
    pub(crate) id: String,
    pub(crate) name: String,
    pub(crate) position: i64,
    pub(crate) reference: String,
    pub(crate) var_type: VariantType,
    pub(crate) alternatives: Option<String>,
    pub(crate) ann: Option<Vec<Vec<String>>>,
    pub(crate) vis: String,
}


pub(crate) fn make_report(vcf_path: &Path, fasta_path: &Path, bam_path: &Path, chrom: String) -> Result<Vec<Report>, Box<dyn Error>> {
    let mut vcf = rust_htslib::bcf::Reader::from_path(&vcf_path).unwrap();
    let header = vcf.header().clone();

    let mut reports = Vec::new();

    for v in vcf.records() {
        let mut variant = v.unwrap();

        let n = header.rid2name(variant.rid().unwrap()).unwrap();
        let i = variant.id();

        let mut name = String::from("");
        for c in n {
            name.push(*c as char);
        }

        let mut id = String::from("");
        for c in i {
            id.push(c as char);
        }

        let pos = variant.pos();
        let end_pos = match variant.info(b"END").integer() {
            Ok(Some(end_pos)) => {
                // Subtraction of 0.5 because of the 0-based positioning in the whole plot
                //TODO: Ask @johanneskoester if the END Info Tag is 0-based or not
                let end_pos = end_pos[0] as f64 - 0.5; // -1 due to 0-basing, + 0.5 du to end pos
                Some(end_pos)
            }
            _ => None,
        };

        let mut ann_strings = Vec::new();
        let ann = variant.info(b"ANN").string()?;

        if let Some(ann) = ann {
            for entry in ann {
                let fields: Vec<_> = entry.split(|c| *c == b'|').collect();
                let mut s = Vec::new();
                for f in fields {
                    let mut attr = String::from("");
                    attr.push_str(std::str::from_utf8(f).unwrap());
                    s.push(attr);
                }
                ann_strings.push(s);
            }
        }


        let alleles = variant.alleles();


        if alleles.len() > 0 {
            let ref_vec = alleles[0].clone();
            let mut rfrce = String::from("");

            let mut len: u8 = 0;

            for c in ref_vec {
                rfrce.push(*c as char);
                len += 1;
            }

            for i in 1..alleles.len() {
                let alt = alleles[i];
                let var_string = String::from("Variant");

                if alt == b"<DEL>" {
                    let var = StaticVariant {
                        marker_type: var_string,
                        reference: rfrce.clone(),
                        alternatives: None,
                        start_position: pos as f64 - 0.5,
                        end_position: end_pos.unwrap(),
                        row: -1,
                        var_type: VariantType::Deletion,
                    };

                    let content = create_report_data(fasta_path, var, bam_path, chrom.clone(), variant.pos() as u64 - 75, end_pos.unwrap() as u64 + 75);
                    let visualization = manipulate_json(content, variant.pos() as u64 - 75, end_pos.unwrap() as u64 + 75);

                    let r = Report {
                        id: id.clone(),
                        name: name.clone(),
                        position: variant.pos(),
                        reference: rfrce.clone(),
                        var_type: VariantType::Deletion,
                        alternatives: None,
                        ann: Some(ann_strings.clone()),
                        vis: visualization.to_string(),
                    };

                    reports.push(r);
                } else if alt == b"<INV>" {
                    let rev: String = rfrce.chars().rev().collect();

                    let var = StaticVariant {
                        marker_type: var_string,
                        reference: rfrce.clone(),
                        alternatives: Some(rev.clone()),
                        start_position: pos as f64 - 0.5,
                        end_position: end_pos.unwrap(),
                        row: -1,
                        var_type: VariantType::Inversion,
                    };

                    let content = create_report_data(fasta_path, var, bam_path, chrom.clone(), variant.pos() as u64 - 75, end_pos.unwrap() as u64 + 75);
                    let visualization = manipulate_json(content, variant.pos() as u64 - 75, end_pos.unwrap() as u64 + 75);

                    let r = Report {
                        id: id.clone(),
                        name: name.clone(),
                        position: variant.pos(),
                        reference: rfrce.clone(),
                        var_type: VariantType::Inversion,
                        alternatives: Some(rev),
                        ann: Some(ann_strings.clone()),
                        vis: visualization.to_string(),
                    };

                    reports.push(r);
                } else if alt == b"<DUP>" {
                    let dup: String = [rfrce.clone(),rfrce.clone()].concat();

                    let var = StaticVariant {
                        marker_type: var_string,
                        reference: rfrce.clone(),
                        alternatives: Some(dup.clone()),
                        start_position: pos as f64 - 0.5,
                        end_position: end_pos.unwrap(),
                        row: -1,
                        var_type: VariantType::Duplicate,
                    };

                    let content = create_report_data(fasta_path, var, bam_path, chrom.clone(), variant.pos() as u64 - 75, end_pos.unwrap() as u64 + 75);
                    let visualization = manipulate_json(content, variant.pos() as u64 - 75, end_pos.unwrap() as u64 + 75);

                    let r = Report {
                        id: id.clone(),
                        name: name.clone(),
                        position: variant.pos(),
                        reference: rfrce.clone(),
                        var_type: VariantType::Duplicate,
                        alternatives: Some(dup.clone()),
                        ann: Some(ann_strings.clone()),
                        vis: visualization.to_string()
                    };

                    reports.push(r);
                } else {
                    let mut allel = String::from("");

                    for c in alt {
                        if *c as char != '<' && *c as char != '>' {
                            allel.push(*c as char);
                        }
                    }

                    if allel.len() == rfrce.len() {
                        let var = StaticVariant {
                            marker_type: var_string,
                            reference: rfrce.clone(),
                            alternatives: Some(allel.clone()),
                            start_position: pos as f64 - 0.5,
                            end_position: pos as f64 - 0.5 + len as f64,
                            row: -1,
                            var_type: VariantType::Variant,
                        };

                        let content = create_report_data(fasta_path, var, bam_path, chrom.clone(), variant.pos() as u64 - 75, variant.pos() as u64 + len as u64 + 75);
                        let visualization = manipulate_json(content, variant.pos() as u64 - 75, variant.pos() as u64 + len as u64 + 75);

                        let r = Report {
                            id: id.clone(),
                            name: name.clone(),
                            position: variant.pos(),
                            reference: rfrce.clone(),
                            var_type: VariantType::Variant,
                            alternatives: Some(allel.clone()),
                            ann: Some(ann_strings.clone()),
                            vis: visualization.to_string()
                        };

                        reports.push(r);
                    } else if allel.len() > rfrce.len() {
                        let var = StaticVariant {
                            marker_type: var_string,
                            reference: rfrce.clone(),
                            alternatives: Some(allel.clone()),
                            start_position: pos as f64, // start end end + 0.5 due to alignment with insertions from bam
                            end_position: pos as f64 + len as f64,
                            row: -1,
                            var_type: VariantType::Insertion,
                        };

                        let content = create_report_data(fasta_path, var, bam_path, chrom.clone(), variant.pos() as u64 - 75, variant.pos() as u64 + len as u64 + 75);
                        let visualization = manipulate_json(content, variant.pos() as u64 - 75, variant.pos() as u64 + len as u64 + 75);

                        let r = Report {
                            id: id.clone(),
                            name: name.clone(),
                            position: variant.pos(),
                            reference: rfrce.clone(),
                            var_type: VariantType::Insertion,
                            alternatives: Some(allel.clone()),
                            ann: Some(ann_strings.clone()),
                            vis: visualization.to_string()
                        };

                        reports.push(r);
                    } else {
                        let var = StaticVariant {
                            marker_type: var_string,
                            reference: rfrce.clone(),
                            alternatives: Some(allel.clone()),
                            start_position: pos as f64 + 0.5, // start position + 1 due to alignment with deletions from bam (example: ref: ACTT alt: A  -> deletion is just CTT)
                            end_position: pos as f64 - 0.5 + len as f64,
                            row: -1,
                            var_type: VariantType::Deletion,
                        };

                        let content = create_report_data(fasta_path, var, bam_path, chrom.clone(), variant.pos() as u64 - 75, variant.pos() as u64 + len as u64 + 75);
                        let visualization = manipulate_json(content, variant.pos() as u64 - 75, variant.pos() as u64 + len as u64 + 75);

                        let r = Report {
                            id: id.clone(),
                            name: name.clone(),
                            position: variant.pos(),
                            reference: rfrce.clone(),
                            var_type: VariantType::Deletion,
                            alternatives: Some(allel.clone()),
                            ann: Some(ann_strings.clone()),
                            vis: visualization.to_string()
                        };

                        reports.push(r);
                    }
                }
            }
        }

    }

    Ok(reports.clone())
}

pub fn create_report_data(fasta_path: &Path, variant: StaticVariant, bam_path: &Path, chrom: String, from: u64, to: u64) -> Json {
    // Daten hier sammeln
    let fasta = json!(read_fasta(fasta_path.clone(), chrom.clone(), from, to));

    let (bases, matches) = get_static_reads(bam_path, fasta_path, chrom.clone(), from, to);

    let alignments1 = json!(bases);
    let alignments2 = json!(matches);


    let mut v = Vec::new();
    v.push(variant);
    let variant = json!(v);

    let fasta_string = fasta.to_string();
    let f = fasta_string.trim_end_matches(']');

    let alignment_string1 = alignments1.to_string();
    let alignment_string2 = alignments2.to_string();

    let mut a1 = alignment_string1.replace('[', "");
    a1 = a1.replace(']',"");

    let mut a2 = alignment_string2.replace('[', "");
    a2 = a2.replace(']',"");

    let variant_string = variant.to_string();
    let v = variant_string.trim_start_matches('[');
    let v_empty = v.trim_end_matches(']');

    static T: &str = ",";

    let r:String;

    if v_empty.is_empty() {
        r = [f,T,&a1,T,&a2,v].concat();
    } else if a1.is_empty() {
        r = [f,T,&a1,&a2,T,v].concat();
    } else {
        r = [f,T,&a1,T,&a2,T,v].concat();
    }

    let values = Json::from_str(&r).unwrap();

    values
}
