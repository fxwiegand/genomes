#![feature(proc_macro_hygiene, decl_macro)]

#[macro_use] extern crate rocket;
#[macro_use] extern crate rust_htslib;
#[macro_use] extern crate bit_vec;
#[macro_use] extern crate rocket_contrib;
//#[macro_use] extern crate rustc_serialize;
#[macro_use] extern crate serde;

#[cfg(test)] mod tests;

use rust_htslib::bam;
use rust_htslib::prelude::*;
use bit_vec::BitVec;
use std::fmt;
//use rustc_serialize::json::{self, ToJson, Json};
use std::collections::BTreeMap;
use rocket_contrib::json::Json;
use serde::ser::{Serialize, SerializeStruct, Serializer};

#[get("/api/v1/alignment/<index>")]
fn genome(index: i32) -> Json<Rd> {Json(read_bam(index))}

fn decode_flags(code :u16) -> String {
    let mut flag_string = String::from("");

    const FLAG_1: &'static str = "1: template having multiple segments in sequencing ";
    const FLAG_2: &'static str = "2: each segment properly aligned according to the aligner ";
    const FLAG_3: &'static str = "4: segment unmapped ";
    const FLAG_4: &'static str = "8: next segment in the template unmapped ";
    const FLAG_5: &'static str = "16: SEQ being reverse complemented ";
    const FLAG_6: &'static str = "32: SEQ of the next segment in the template being reverse complemented ";
    const FLAG_7: &'static str = "64: the first segment in the template ";
    const FLAG_8: &'static str = "128: the last segment in the template ";
    const FLAG_9: &'static str = "256: secondary alignment ";
    const FLAG_10: &'static str = "512: not passing filters, such as platform/vendor quality controls ";
    const FLAG_11: &'static str = "1024:PCR or optical duplicate ";
    const FLAG_12: &'static str = "2048: supplementary alignment ";

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

fn read_bam (index: i32) -> Rd {
    let mut bam = bam::Reader::from_path(&"data/test.bam").unwrap();
    let header = bam::Header::from_template(bam.header());

    let mut read:Rd = Rd {
        header: String::from(""),
        sequenz: String::from(""),
        position: 0,
        cigar_str: String::from(""),
        flags: String::from(""),
        name: String::from(""),
    };

    for r in bam.records() {
        let record = r.unwrap();
        let head = header.to_bytes();

        //Cigar String
        let cigstring = record.cigar().to_string();

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

        //Header
        let mut hd = String::from("");
        for b in head {
            hd.push(b as char);
        }

        //Name
        let n = record.qname();
        let mut name = String::from("");
        for a in n {
            name.push(*a as char);
        }

        let read1 = Rd {
            header: hd,
            sequenz: sequenz,
            position: pos,
            cigar_str: cigstring,
            flags: flag_string,
            name: name,
        };

        if index == read1.position {
            read = read1;
        }
    }
    read
}

fn main() {
    rocket::ignite().mount("/", routes![genome]).launch();
}

#[derive(Serialize)]
struct Rd {
    header: String,
    sequenz: String,
    position: i32,
    cigar_str: String,
    flags: String,
    name: String,
}

impl fmt::Display for Rd {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {}, {}, {}, {}, {})", self.header, self.sequenz, self.position, self.cigar_str, self.flags,
        self.name)
    }
}
