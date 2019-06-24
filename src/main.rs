#![feature(proc_macro_hygiene, decl_macro)]

#[macro_use] extern crate rocket;
#[macro_use] extern crate serde;

extern crate rocket_contrib;
extern crate bit_vec;
extern crate bio;
extern crate rust_htslib;

mod alignment_reader;
mod fasta_reader;
mod variant_reader;

#[cfg(test)] mod tests;

use std::path::Path;
use std::env;
use rocket_contrib::json::Json;
use rocket_contrib::serve::StaticFiles;
use rocket::State;
use alignment_reader::{get_reads, AlignmentNucleobase, read_indexed_bam};
use alignment_reader::Alignment;
use fasta_reader::Nucleobase;
use fasta_reader::read_fasta;
use variant_reader::read_vcf;
use variant_reader::read_indexed_vcf;
use variant_reader::Variant;


#[get("/reference/<chromosome>/<from>/<to>")]
fn reference(args: State<Vec<String>>, chromosome: u8, from: u64, to: u64) -> Json<Vec<Nucleobase>> {Json(read_fasta(Path::new(&args[2].clone()), chromosome - 1, from, to))}

#[get("/alignment/<chromosome>/<from>/<to>")]
fn alignment(args: State<Vec<String>>, chromosome: u8, from: u32, to: u32) -> Json<Vec<AlignmentNucleobase>> {
    Json(get_reads(Path::new(&args[1].clone()), chromosome, from, to))
}
#[get("/variant/<chromosome>/<from>/<to>")]
fn variant(args: State<Vec<String>>, chromosome: u8, from: u32, to: u32) -> Json<Vec<Variant>> {
    Json(read_indexed_vcf(Path::new(&args[3].clone()), chromosome, from, to))
}

fn main() {

    let args: Vec<String> = env::args().collect();

    let bam_filename = args[1].clone();
    let bam_path = Path::new(&bam_filename);

    let vcf_filename = args[3].clone();
    let vcf_path = Path::new(&vcf_filename);



    rocket::ignite()
        .manage(args)
        .mount("/",  StaticFiles::from("client"))
        .mount("/api/v1", routes![reference, alignment, variant])
        .launch();

}


