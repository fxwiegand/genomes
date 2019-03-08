#![feature(proc_macro_hygiene, decl_macro)]

#[macro_use] extern crate rocket;
#[macro_use] extern crate rust_htslib;


#[cfg(test)] mod tests;

use rust_htslib::bam;
use rust_htslib::prelude::*;

#[get("/lisbert/<verb>")]
fn hello(verb: String) -> String {format!("Lisa {}!", verb)}

fn main() {
    rocket::ignite().mount("/", routes![hello]).launch();

    let mut bam = bam::Reader::from_path(&"data/test.bam").unwrap();
    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_path(&"data/out.bam", &header).unwrap();

    for r in bam.records() {
        let record = r.unwrap();
        if record.is_reverse() {
            out.write(&record).unwrap();
        }
    }
}


