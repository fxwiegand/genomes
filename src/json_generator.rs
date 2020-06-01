use std::path::Path;
use rustc_serialize::json::Json;
use fasta_reader::read_fasta;
use static_reader::{get_static_reads,get_static_variants};
use serde_json::{Value};
use std::fs;


pub fn create_data(fasta_path: &Path, vcf_path: &Path, bam_path: &Path, chrom: String, from: u64, to: u64) -> Json {
    let mut data = Vec::new();

    for f in read_fasta(fasta_path.clone(), chrom.clone(), from, to) {
        let nucleobase = json!(f);
        data.push(nucleobase);
    }
    let (bases, matches) = get_static_reads(bam_path, fasta_path, chrom.clone(), from, to);

    for b in bases {
        let base = json!(b);
        data.push(base);
    }

    for m in matches {
        let mat = json!(m);
        data.push(mat);
    }

    for v in get_static_variants(vcf_path, chrom.clone(), from, to) {
        let variant = json!(v);
        data.push(variant);
    }

    let values = Json::from_str(&json!(data).to_string()).unwrap();

    values
}



pub fn manipulate_json(data: Json, from: u64, to: u64) -> Value {
    let vega_path = Path::new("./static/vegaSpecs.json");
    let json_string = fs::read_to_string(vega_path).unwrap();

    let mut vega_specs :Value = serde_json::from_str(&json_string).unwrap();
    let values: Value = serde_json::from_str(&data.to_string()).unwrap();
    let mut values = json!({"values": values, "name": "fasta"});

    let v = values["values"].as_array().unwrap().clone();


    for i in 0..v.len() {
        let k = v[i]["marker_type"].clone().as_str().unwrap().to_owned();

        if k == "A" || k == "T" || k == "G" || k == "C" || k == "U" {
            values["values"][i]["base"] = values["values"][i]["marker_type"].clone();
        } else if k == "Deletion" || k == "Match" || k == "Pairing" || k == "Duplicate" || k == "Inversion" {
            values["values"][i]["typ"] = values["values"][i]["marker_type"].clone();
        } else if k == "Insertion"{
            values["values"][i]["typ"] = values["values"][i]["marker_type"].clone();
            values["values"][i]["inserts"] = values["values"][i]["bases"].clone();
        }
    }

    vega_specs["width"] = json!(700);
    let domain = json!([from, to]);

    vega_specs["scales"][0]["domain"] = domain;
    vega_specs["data"][1] = values;

    vega_specs
}