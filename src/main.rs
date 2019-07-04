#![feature(proc_macro_hygiene, decl_macro)]

#[macro_use] extern crate rocket;
#[macro_use] extern crate serde;
#[macro_use] extern crate serde_json;

extern crate rocket_contrib;
extern crate bit_vec;
extern crate bio;
extern crate rust_htslib;
extern crate jsonm;

mod alignment_reader;
mod fasta_reader;
mod variant_reader;

#[cfg(test)] mod tests;

use std::path::Path;
use std::env;
use rocket_contrib::json::{Json};
use rocket_contrib::serve::StaticFiles;
use rocket::State;
use jsonm::packer::{PackOptions, Packer};
use alignment_reader::{get_reads};
use fasta_reader::read_fasta;
use variant_reader::read_indexed_vcf;
use serde_json::Value;


#[get("/reference/<chromosome>/<from>/<to>")]
fn reference(args: State<Vec<String>>, chromosome: String, from: u64, to: u64) -> Json<Value> {
    let response = read_fasta(Path::new(&args[2].clone()), chromosome, from, to);
    let mut packer = Packer::new();
    let options = PackOptions::new();
    let packed = packer.pack(&json!(response), &options).unwrap();
    Json(packed)
}

#[get("/alignment/<chromosome>/<from>/<to>")]
fn alignment(args: State<Vec<String>>, chromosome: String, from: u32, to: u32) -> Json<Value> {
    let response = get_reads(Path::new(&args[1].clone()), Path::new(&args[2].clone()) , chromosome, from, to);
    let mut packer = Packer::new();
    let options = PackOptions::new();
    let packed = packer.pack(&json!(response), &options).unwrap();
    Json(packed)
}

#[get("/variant/<chromosome>/<from>/<to>")]
fn variant(args: State<Vec<String>>, chromosome: String, from: u32, to: u32) -> Json<Value> {
    let response = read_indexed_vcf(Path::new(&args[3].clone()), chromosome, from, to);
    let mut packer = Packer::new();
    let options = PackOptions::new();
    let packed = packer.pack(&json!(response), &options).unwrap();
    Json(packed)
}

fn main() {

    let args: Vec<String> = env::args().collect();

    rocket::ignite()
        .manage(args)
        .mount("/",  StaticFiles::from("client"))
        .mount("/api/v1", routes![reference, alignment, variant])
        .launch();

}


