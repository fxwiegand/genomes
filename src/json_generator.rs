use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use rustc_serialize::json::Json;
use fasta_reader::read_fasta;
use variant_reader::read_indexed_vcf;
use alignment_reader::get_static_reads;


pub fn create_data(fasta_path: &Path, vcf_path: &Path, bam_path: &Path, chrom: String, from: u32, to: u32) -> std::io::Result<()> {
    // Daten hier sammeln
    let fasta = json!(read_fasta(fasta_path.clone(), chrom.clone(), from as u64, to as u64));
    let alignments = json!(get_static_reads(bam_path, fasta_path, chrom.clone(), from, to));
    let variant = json!(read_indexed_vcf(vcf_path, chrom.clone(), from, to));


    let fasta_string = fasta.to_string();
    let f = fasta_string.trim_end_matches(']');


    let alignment_string = alignments.to_string();
    let mut a = alignment_string.trim_start_matches('[');
    a = a.trim_end_matches(']');


    let variant_string = variant.to_string();
    let v = variant_string.trim_start_matches('[');
    let v_empty = v.trim_end_matches(']');

    static T: &str = ",";

    let mut r= String::from("[]");

    if v_empty.is_empty() {
        r = [f,T,a,v].concat();
    } else if a.is_empty() {
        r = [f,T,v].concat();
    } else {
        r = [f,T,a,T,v].concat();
    }

    let values = Json::from_str(&r).unwrap();


    let mut file = File::create("data.json")?;
    file.write_all(values.to_string().as_bytes())?;
    Ok(())
}