#![feature(proc_macro_hygiene, decl_macro)]

#[macro_use] extern crate serde;
extern crate rust_htslib;
extern crate bit_vec;

use rust_htslib::{bam};
use rust_htslib::prelude::*;
use std::fmt;
use std::path::Path;
use std::env;
use std::collections::BTreeMap;

#[derive(Serialize, Clone)]
pub struct Alignment {
    index: i32,
    header: String,
    sequence: String,
    pos: i32,
    cigar: String,
    flags: BTreeMap<u16, &'static str>,
    name: String,
}

impl fmt::Display for Alignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut flag_string = String::from("");
        for (_key, flag) in &self.flags {
            flag_string.push_str(flag);
        }
        write!(f, "({}, {}, {}, {}, {}, {})", self.header, self.sequence, self.pos, self.cigar, flag_string,
               self.name)
    }
}


pub fn decode_flags(code :u16) -> BTreeMap<u16, &'static str> {
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

    let mut flags_map = BTreeMap::new();
    flags_map.insert(0x1, FLAG_1);
    flags_map.insert(0x2, FLAG_2);
    flags_map.insert(0x4, FLAG_3);
    flags_map.insert(0x8, FLAG_4);
    flags_map.insert(0x10, FLAG_5);
    flags_map.insert(0x20, FLAG_6);
    flags_map.insert(0x40, FLAG_7);
    flags_map.insert(0x80, FLAG_8);
    flags_map.insert(0x100, FLAG_9);
    flags_map.insert(0x200, FLAG_10);
    flags_map.insert(0x400, FLAG_11);
    flags_map.insert(0x800, FLAG_12);

    for (flag, text) in flags_map {
        if (flag & code) == flag {
            string_map.insert(flag, text);
        }
    }

    string_map
}

pub fn count_alignments(path: &Path)-> u32 {
    let mut bam = bam::Reader::from_path(path).unwrap();
    //let header = bam::Header::from_template(bam.header());
    let mut count:u32= 0;
    for _r in bam.records() {
        count += 1;
    }

    count
}


pub fn read_bam(path: &Path) -> Vec<Alignment> {
    let mut bam = bam::Reader::from_path(path).unwrap();
    let header = bam::Header::from_template(bam.header());

    let mut alignments:Vec<Alignment> = Vec::new();
    let mut ind = 0;

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

        let read = Alignment {
            index: ind,
            header: hd,
            sequence: sequenz,
            pos: pos,
            cigar: cigstring,
            flags: flag_string,
            name: name,
        };

        alignments.push(read);
        ind +=1;
    }

    alignments

}