#![feature(proc_macro_hygiene, decl_macro)]

#[macro_use] extern crate rocket;
#[macro_use] extern crate serde;
#[macro_use] extern crate serde_json;
#[macro_use] extern crate log;

extern crate rocket_contrib;
extern crate bit_vec;
extern crate bio;
extern crate rust_htslib;
extern crate rustc_serialize;
extern crate regex;
extern crate clap;
extern crate tera;

mod alignment_reader;
mod fasta_reader;
mod variant_reader;
mod json_generator;
mod static_reader;

#[cfg(test)] mod variant_tests;
#[cfg(test)] mod reference_tests;
#[cfg(test)] mod alignment_tests;

use std::path::Path;
use std::process::Command;
use std::str::FromStr;
use std::io::{self, Write};
use rocket_contrib::json::{Json};
use rocket_contrib::serve::StaticFiles;
use rocket_contrib::compression::Compression;
use rocket::State;
use clap::{Arg, App, SubCommand, ArgMatches};
use alignment_reader::{get_reads, AlignmentNucleobase, AlignmentMatch};
use fasta_reader::{read_fasta, Nucleobase};
use variant_reader::{read_indexed_vcf, read_vcf,Variant};
use json_generator::create_data;
use rocket_contrib::templates::Template;
use tera::Context;
use std::collections::HashMap;


#[get("/reference/<chromosome>/<from>/<to>")]
fn reference(params: State<ArgMatches>, chromosome: String, from: u64, to: u64) -> Json<Vec<Nucleobase>> {
    let response = read_fasta(Path::new(params.value_of("fasta file").unwrap()), chromosome, from, to);
    Json(response)
}

#[get("/alignment/<chromosome>/<from>/<to>")]
fn alignment(params: State<ArgMatches>, chromosome: String, from: u64, to: u64) -> Json<(Vec<AlignmentNucleobase>,Vec<AlignmentMatch>)> {
    let response = get_reads(Path::new(params.value_of("bam file").unwrap()), Path::new(params.value_of("fasta file").unwrap()) , chromosome, from, to);
    Json(response)
}

#[get("/variant/<chromosome>/<from>/<to>")]
fn variant(params: State<ArgMatches>, chromosome: String, from: u64, to: u64) -> Json<Vec<Variant>> {
    let response = read_indexed_vcf(Path::new(params.value_of("vcf file").unwrap()), chromosome, from, to);
    Json(response)
}

#[get("/")]
fn index(params: State<ArgMatches>) -> Template {
    let mut context = HashMap::new();
    context.insert("variants", read_vcf(Path::new(params.value_of("vcf file").unwrap())));

    Template::render("report", &context)
}

fn main() {
    let matches = App::new("gensbock")
        .version("1.0")
        .author("Felix W. <fxwiegand@wgdnet.de>")
        .about("genome viewing in rust")
        .subcommand(SubCommand::with_name("server")
            .about("starts server")
            .version("1.0")
            .author("Felix W. <fxwiegand@wgdnet.de>")
            .arg(Arg::with_name("bam file")
                .required(true)
                .help("your input bam file")
                .index(1))
            .arg(Arg::with_name("fasta file")
                .required(true)
                .help("your input fasta file")
                .index(2))
            .arg(Arg::with_name("vcf file")
                .required(true)
                .help("your input vcf file")
                .index(3)))
        .subcommand(SubCommand::with_name("static")
            .about("outputs vega specs")
            .version("1.0")
            .author("Felix W. <fxwiegand@wgdnet.de>")
            .arg(Arg::with_name("bam file")
                .required(true)
                .help("your input bam file")
                .index(1))
            .arg(Arg::with_name("fasta file")
                .required(true)
                .help("your input fasta file")
                .index(2))
            .arg(Arg::with_name("vcf file")
                .required(true)
                .help("your input vcf file")
                .index(3))
            .arg(Arg::with_name("chromosome")
                .required(true)
                .help("the chromosome you want to visualize")
                .index(4))
            .arg(Arg::with_name("from")
                .required(true)
                .help("the start of the region you want to visualize")
                .index(5))
            .arg(Arg::with_name("to")
                .required(true)
                .help("the end of the region you want to visualize")
                .index(6)))
        .subcommand(SubCommand::with_name("report")
        .arg(Arg::with_name("vcf file")
            .required(true)
            .help("your input vcf file")
            .index(1)))
        .get_matches();

    match matches.subcommand_name() {
        Some("server") => {
            let params = matches.subcommand_matches("server").unwrap().clone();

            rocket::ignite()
                .manage(params)
                .mount("/",  StaticFiles::from("client"))
                .mount("/api/v1", routes![reference, alignment, variant])
                .attach(Compression::fairing())
                .launch();
        },
        Some("static") => {
            let static_matches = matches.subcommand_matches("static").unwrap();

            let fasta_path = Path::new(static_matches.value_of("fasta file").unwrap());
            let bam_path = Path::new(static_matches.value_of("bam file").unwrap());
            let vcf_path = Path::new(static_matches.value_of("vcf file").unwrap());
            let chromosome = String::from(static_matches.value_of("chromosome").unwrap());
            let from = u64::from_str(static_matches.value_of("from").unwrap()).unwrap();
            let to = u64::from_str(static_matches.value_of("to").unwrap()).unwrap();


            let _data = create_data(&fasta_path,&vcf_path,&bam_path,chromosome,from,to);
            let a :String = from.to_string();
            let b :String = to.to_string();
            let py = String::from("src/jsonMerge.py");

            let output = {
                Command::new("python")
                    .arg(py)
                    .arg(a)
                    .arg(b)
                    .output()
                    .expect("failed to execute process")
            };

            let msg = output.stdout;
            io::stdout().write(msg.as_ref()).unwrap();
        },
        Some("report") => {
            let params = matches.subcommand_matches("report").unwrap().clone();

            rocket::ignite()
                .manage(params)
                .mount("/", routes![index])
                .attach(Template::fairing())
                .launch();
        },
        None        => println!("Try using a subcommand. Type help for more."),
        _           => unreachable!(), // Assuming you've listed all direct children above, this is unreachable
    }
}


