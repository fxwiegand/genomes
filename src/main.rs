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
fn genome(index: i32) -> Json<Alignment> {Json(read_bam(index))}

#[get("/api/v1/count")]
fn count() -> Json<u32> {Json(count_alignments())}

fn decode_flags(code :u16) -> BTreeMap<u32, &'static str> {
    let mut string_map = BTreeMap::new();

    const FLAG_1: &'static str = "template having multiple segments in sequencing";
    const FLAG_2: &'static str = "each segment properly aligned according to the aligner";
    const FLAG_3: &'static str = "segment unmapped";
    const FLAG_4: &'static str = "next segment in the template unmapped";
    const FLAG_5: &'static str = "SEQ being reverse complemented";
    const FLAG_6: &'static str = "SEQ of the next segment in the template being reverse complemented";
    const FLAG_7: &'static str = "the first segment in the template ";
    const FLAG_8: &'static str = "the last segment in the template";
    const FLAG_9: &'static str = "secondary alignment";
    const FLAG_10: &'static str = "not passing filters, such as platform/vendor quality controls";
    const FLAG_11: &'static str = "PCR or optical duplicate";
    const FLAG_12: &'static str = "supplementary alignment";

    if (0x1 & code) == 0x1 {
        string_map.insert(1 as u32, FLAG_1);
    }
    if  (0x2 & code) == 0x2  {
        string_map.insert(2 as u32, FLAG_2);
    }
    if (0x4 & code) == 0x4 {
        string_map.insert(4 as u32, FLAG_3);
    }
    if  (0x8 & code) == 0x8  {
        string_map.insert(8 as u32, FLAG_4);
    }
    if (0x10 & code) == 0x10 {
        string_map.insert(16 as u32, FLAG_5);
    }
    if  (0x20 & code) == 0x20  {
        string_map.insert(32 as u32, FLAG_6);
    }
    if (0x40 & code) == 0x40 {
        string_map.insert(64 as u32, FLAG_7);
    }
    if  (0x80 & code) == 0x80  {
        string_map.insert(128 as u32, FLAG_8);
    }
    if  (0x100 & code) == 0x100  {
        string_map.insert(256 as u32, FLAG_9);
    }
    if (0x200 & code) == 0x200 {
        string_map.insert(512 as u32, FLAG_10);
    }
    if  (0x400 & code) == 0x400  {
        string_map.insert(1024 as u32, FLAG_11);
    }
    if (0x800 & code) == 0x800 {
        string_map.insert(2048 as u32, FLAG_12);
    }

    string_map
}

fn count_alignments()-> u32 {
    let mut bam = bam::Reader::from_path(&"data/test.bam").unwrap();
    let header = bam::Header::from_template(bam.header());
    let mut count:u32= 0;
    for r in bam.records() {
        count += 1;
    }

    count
}

fn read_bam(index: i32) -> Alignment {
    let mut bam = bam::Reader::from_path(&"data/test.bam").unwrap();
    let header = bam::Header::from_template(bam.header());

    let mut read:Alignment = Alignment {
        header: String::from(""),
        sequenz: String::from(""),
        position: 0,
        cigar_str: String::from(""),
        flags: BTreeMap::new(),
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

        let read1 = Alignment {
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
    rocket::ignite().mount("/", routes![genome, count]).launch();
}

#[derive(Serialize)]
struct Alignment {
    header: String,
    sequenz: String,
    position: i32,
    cigar_str: String,
    flags: BTreeMap<u32, &'static str>,
    name: String,
}

impl fmt::Display for Alignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut flag_string = String::from("");
        for (_key, flag) in &self.flags {
            flag_string.push_str(flag);
        }
        write!(f, "({}, {}, {}, {}, {}, {})", self.header, self.sequenz, self.position, self.cigar_str, flag_string,
        self.name)
    }
}
