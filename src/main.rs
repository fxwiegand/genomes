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
use fasta_reader::Nucleobase;
use fasta_reader::read_fasta;
use std::collections::BTreeMap;


#[get("/alignment/<index>")]
fn genome(index: usize, alignments: State<Vec<Alignment>>) -> Json<Alignment> {Json(alignments[index].clone())}

#[get("/count")]
fn count(size: State<u32>) -> Json<u32> {Json(size.clone())}

#[get("/reference/<chromosome>/<from>/<to>")]
fn reference(args: State<Vec<String>>, chromosome: u8, from: u64, to: u64) -> Json<Vec<Nucleobase>> {Json(read_fasta(Path::new(&args[2].clone()), chromosome - 1, from, to))}

fn main() {
    let args: Vec<String> = env::args().collect();

    // The first argument is the path that was used to call the program.
    ////println!("I got {:?} arguments: {:?}.", args.len() - 1, &args[1..]);


    let bam_filename = args[1].clone();
    //let fasta_filename   = args[2].clone();
    let bam_path = Path::new(&bam_filename);
    //let fasta_path: &Path = Path::new(&fasta_filename);



    let alignments = read_bam(bam_path);

    let size = count_alignments(bam_path);

    let _index = Route::ranked(1, Method::Get, "/", StaticFiles::from("client"));



    rocket::ignite()
        .manage(alignments)
        .manage(size)
        .manage(args)
        .mount("/home",  StaticFiles::from("client"))
        .mount("/static", StaticFiles::from("client"))
        .mount("/optics", StaticFiles::from("optics"))
        .mount("/api/v1", routes![genome, count, reference])
        .launch();

}


