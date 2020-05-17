use std::path::Path;
use rustc_serialize::json::Json;
use fasta_reader::read_fasta;
use static_reader::{get_static_reads,get_static_variants};
use serde_json::{Value};
use std::fs;


pub fn create_data(fasta_path: &Path, vcf_path: &Path, bam_path: &Path, chrom: String, from: u64, to: u64) -> Json {
    // Daten hier sammeln
    let fasta = json!(read_fasta(fasta_path.clone(), chrom.clone(), from, to));
    let alignments = json!(get_static_reads(bam_path, fasta_path, chrom.clone(), from, to));
    let variant = json!(get_static_variants(vcf_path, chrom.clone(), from, to));

    let fasta_string = fasta.to_string();
    let f = fasta_string.trim_end_matches(']');

    let alignment_string = alignments.to_string();
    let mut a = alignment_string.replace('[', "");
    a = a.replace(']',"");


    let variant_string = variant.to_string();
    let v = variant_string.trim_start_matches('[');
    let v_empty = v.trim_end_matches(']');

    static T: &str = ",";

    let r:String;

    if v_empty.is_empty() {
        r = [f,T,&a,v].concat();
    } else if a.is_empty() {
        r = [f,T,v].concat();
    } else {
        r = [f,T,&a,T,v].concat();
    }

    let values = Json::from_str(&r).unwrap();

    values
}



pub fn manipulate_json(data: Json, from: u64, to: u64) -> Value {
    let vega_path = Path::new("./client/vegaSpecs.json");
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