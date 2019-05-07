#![feature(proc_macro_hygiene, decl_macro)]

#[macro_use] extern crate rocket;
#[macro_use] extern crate serde;

extern crate rocket_contrib;
extern crate rust_htslib;
extern crate bit_vec;
extern crate bio;

mod alignment_reader;
mod fasta_reader;

#[cfg(test)] mod tests;

use std::path::Path;
use std::env;
use rocket_contrib::json::Json;
use rocket_contrib::serve::StaticFiles;
use rocket::{State, Route};
use rocket::http::Method;
use alignment_reader::count_alignments;
use alignment_reader:: read_bam;
use alignment_reader::Alignment;
use fasta_reader::read_fasta;




#[get("/alignment/<index>")]
fn genome(index: usize, alignments: State<Vec<Alignment>>) -> Json<Alignment> {Json(alignments[index].clone())}

#[get("/count")]
fn count(size: State<u32>) -> Json<u32> {Json(size.clone())}

fn main() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let path = Path::new(filename);

    let alignments = read_bam(path);

    let size = count_alignments(path);


    let _index = Route::ranked(1, Method::Get, "/", StaticFiles::from("client"));

    let fasta_path = Path::new("data/human_b36_male.fa");

    let fasta = read_fasta(fasta_path, 0, 1, 100);

    println!("{}", fasta);

//    rocket::ignite()
//        .manage(alignments)
//        .manage(size)
//        .mount("/home",  StaticFiles::from("client"))
//        .mount("/static", StaticFiles::from("client"))
//        .mount("/optics", StaticFiles::from("optics"))
//        .mount("/api/v1", routes![genome, count])
//        .launch();

}


