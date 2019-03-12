#![feature(proc_macro_hygiene, decl_macro)]

#[macro_use] extern crate rocket;
#[macro_use] extern crate rust_htslib;
#[macro_use] extern crate bit_vec;

#[cfg(test)] mod tests;

use rust_htslib::bam;
use rust_htslib::prelude::*;
use bit_vec::BitVec;


#[get("/lisbert/<verb>")]
fn hello(verb: String) -> String {format!("Lisa {}!", verb)}

fn decode_flags(code :u16) -> String {
    let mut flag_string = String::from("");

    let FLAG_1 = String::from("1: template having multiple segments in sequencing ");
    let FLAG_2 = String::from("2: each segment properly aligned according to the aligner ");
    let FLAG_3 = String::from("4: segment unmapped ");
    let FLAG_4 = String::from("8: next segment in the template unmapped ");
    let FLAG_5 = String::from("16: SEQ being reverse complemented ");
    let FLAG_6 = String::from("32: SEQ of the next segment in the template being reverse complemented ");
    let FLAG_7 = String::from("64: the first segment in the template ");
    let FLAG_8 = String::from("128: the last segment in the template ");
    let FLAG_9 = String::from("256: secondary alignment ");
    let FLAG_10 = String::from("512: not passing filters, such as platform/vendor quality controls ");
    let FLAG_11 = String::from("1024:PCR or optical duplicate ");
    let FLAG_12 = String::from("2048: supplementary alignment ");

    if (0x1 & code) == 0x1 {
        flag_string.push_str(&FLAG_1);
    }
    if  (0x2 & code) == 0x2  {
        flag_string.push_str(&FLAG_2);
    }
    if (0x4 & code) == 0x4 {
        flag_string.push_str(&FLAG_3);
    }
    if  (0x8 & code) == 0x8  {
        flag_string.push_str(&FLAG_4);
    }
    if (0x10 & code) == 0x10 {
        flag_string.push_str(&FLAG_5);
    }
    if  (0x20 & code) == 0x20  {
        flag_string.push_str(&FLAG_6);
    }
    if (0x40 & code) == 0x40 {
        flag_string.push_str(&FLAG_7);
    }
    if  (0x80 & code) == 0x80  {
        flag_string.push_str(&FLAG_8);
    }
    if  (0x100 & code) == 0x100  {
        flag_string.push_str(&FLAG_9);
    }
    if (0x200 & code) == 0x200 {
        flag_string.push_str(&FLAG_10);
    }
    if  (0x400 & code) == 0x400  {
        flag_string.push_str(&FLAG_11);
    }
    if (0x800 & code) == 0x800 {
        flag_string.push_str(&FLAG_12);
    }

    flag_string
}


fn main() {
    //rocket::ignite().mount("/", routes![hello]).launch();

    let mut bam = bam::Reader::from_path(&"data/test.bam").unwrap();
    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_path(&"data/out.bam", &header).unwrap();

    for r in bam.records() {
        let record = r.unwrap();
        let head = header.to_bytes();

        //Cigar String
        let cigstring = record.cigar();

        //Position
        let pos = record.pos();

        //Sequenz
        let seq = record.seq().as_bytes();
        let mut sequenz = String::from("");
        for b in seq {
            sequenz.push(b as char);
        }

        //Flags
        let flgs = record.flags();
        let flag_string = decode_flags(flgs);
        //let bytes = format!()
        //let bv = BitVec::from_bytes();

        //Header
        let mut hd = String::from("");
        for b in head {
            hd.push(b as char);
        }

        println!("Header: {}", hd);
        println!("Sequenz: {}", sequenz);
        println!("Position: {}", pos);
        println!("CIGAR-String: {}", cigstring);
        println!("Flags: {}", flag_string)
        //println!("Was bin ich?: {}", name as usize);
    }
}


