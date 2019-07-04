extern crate rust_htslib;
extern crate bit_vec;

use rust_htslib::{bam};
use rust_htslib::prelude::*;
use std::fmt;
use std::path::Path;
use std::collections::BTreeMap;
use rust_htslib::bam::record::CigarStringView;
use fasta_reader::{read_fasta};

#[derive(Serialize, Clone)]
pub enum Marker {
    A,
    T,
    C,
    G,
    N,
    Deletion,
    Insertion,
    Match,
}

#[derive(Clone)]
pub struct Alignment {
    sequence: String,
    pos: i32,
    length: u16,
    flags: BTreeMap<u16, &'static str>,
    name: String,
    cigar: CigarStringView,
    paired: bool,
    mate_pos: i32,
}

#[derive(Serialize, Clone)]
pub struct AlignmentNucleobase {
    marker_type: Marker,
    bases: String,
    position: f32,
    flags: BTreeMap<u16, &'static str>,
    name: String,
    read_start: u32,
    read_end: u32,
}


impl fmt::Display for Alignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut flag_string = String::from("");
        for (_key, flag) in &self.flags {
            flag_string.push_str(flag);
        }
        write!(f, "({}, {}, {}, {}, {})", self.sequence, self.pos, self.cigar, flag_string,
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
    const FLAG_12: &'static str = "vega lite lines";

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

pub fn read_indexed_bam(path: &Path, chrom: String, from: u32, to: u32) -> Vec<Alignment> {
    let mut bam = bam::IndexedReader::from_path(&path).unwrap();
    let tid = bam.header().tid(chrom.as_bytes()).unwrap();

    let mut alignments: Vec<Alignment> = Vec::new();

    bam.fetch(tid, from, to).unwrap();

    for r in bam.records() {

        let rec = r.unwrap();

        let a = make_alignment(rec);

        alignments.push(a);
    }

    alignments
}


pub fn read_bam(path: &Path) -> Vec<Alignment> {
    let mut bam = bam::Reader::from_path(path).unwrap();
    let header = bam::Header::from_template(bam.header());

    let mut alignments:Vec<Alignment> = Vec::new();

    for r in bam.records() {
        let record = r.unwrap();
        let _head = header.to_bytes();

        let read = make_alignment(record);

        alignments.push(read);
    }

    alignments

}


fn make_alignment(record: bam::Record) -> Alignment {

    let has_pair = record.is_paired();

    let mate_pos = record.mpos();

    //Cigar String
    let cigstring = record.cigar();

    //Position
    let pos = record.pos();

    //Länge
    let le = record.seq().len() as u16;

    //Sequenz
    let seq = record.seq().as_bytes();
    let mut sequenz = String::from("");
    for b in seq {
        sequenz.push(b as char);
    }

    //Flags
    let flgs = record.flags();
    let flag_string = decode_flags(flgs);

    //Name
    let n = record.qname();
    let mut name = String::from("");
    for a in n {
        name.push(*a as char);
    }

    let read = Alignment {
        sequence: sequenz,
        pos: pos,
        length: le,
        cigar: cigstring,
        flags: flag_string,
        name: name,
        paired: has_pair,
        mate_pos: mate_pos,
    };

    read
}

fn make_nucleobases(fasta_path: &Path, chrom: String, snippets: Vec<Alignment>, from: u32, to: u32) -> Vec<AlignmentNucleobase> {
    let mut bases: Vec<AlignmentNucleobase> = Vec::new();

    let ref_bases = read_fasta(fasta_path, chrom, from as u64, to as u64);

    for s in snippets {
        let mut cigar_offset: i32 = 0;
        let mut read_offset: i32 = 0;
        let base_string = s.sequence.clone();
        let char_vec: Vec<char> = base_string.chars().collect();

        let mut soft_clip_begin = true;

        for c in s.cigar.iter() {
            match c {
                rust_htslib::bam::record::Cigar::Match(c) => {
                    for _i in 0..rust_htslib::bam::record::Cigar::Match(*c).len() {
                        let snip = s.clone();
                        let b = char_vec[cigar_offset as usize];

                        if snip.pos + read_offset >= from as i32 && snip.pos + read_offset < to as i32 {
                            let ref_index = snip.pos + read_offset - from as i32;
                            let ref_base = &ref_bases[ref_index as usize];


                            let m: Marker;
                            if ref_base.get_marker_type() == b {
                                m = Marker::Match; // Match with reference fasta
                            } else {
                                match b { // Mismatch
                                    'A' => m = Marker::A,
                                    'T' => m = Marker::T,
                                    'C' => m = Marker::C,
                                    'N' => m = Marker::N,
                                    'G' => m = Marker::G,
                                    _ => m = Marker::Deletion,
                                }
                            }


                            let p = snip.pos as i32 + read_offset;
                            let f = snip.flags;
                            let n = snip.name;

                            let rs: i32;
                            let re: i32;

                            if snip.paired {
                                if snip.pos < snip.mate_pos {
                                    re = snip.mate_pos + 100;
                                    rs = snip.pos;
                                } else {
                                    rs = snip.mate_pos;
                                    re = snip.pos as i32 + snip.length as i32;
                                }
                            } else {
                                rs = snip.pos;
                                re = snip.pos as i32 + snip.length as i32;
                            }


                            let base = AlignmentNucleobase {
                                marker_type: m,
                                bases: b.to_string(),
                                position: p as f32,
                                flags: f,
                                name: n,
                                read_start: rs as u32,
                                read_end: re as u32,
                            };


                            bases.push(base);
                        }
                        cigar_offset += 1;
                        read_offset += 1;

                    }

                    soft_clip_begin = false;

                }
                rust_htslib::bam::record::Cigar::Ins(c) => {
                    let snip = s.clone();
                    let p: f32 = snip.pos as f32 + read_offset as f32 - 0.5;
                    let m: Marker = Marker::Insertion;
                    let rs = snip.pos;
                    let re = snip.pos as i32 + snip.length as i32;


                    let mut b = String::from("");
                    for i in 0..rust_htslib::bam::record::Cigar::Ins(*c).len() {

                        let char = char_vec[cigar_offset as usize + i as usize];
                        b.push(char);

                    }

                    cigar_offset += 1;

                    let base = AlignmentNucleobase {
                        marker_type: m,
                        bases: b,
                        position: p,
                        flags: snip.flags,
                        name: snip.name,
                        read_start: rs as u32,
                        read_end: re as u32,
                    };



                    if from as f32 <= base.position && base.position <= to as f32 {
                        bases.push(base);
                    }

                    soft_clip_begin = false;

                }
                rust_htslib::bam::record::Cigar::Del(c) => {
                    for _i in 0..rust_htslib::bam::record::Cigar::Del(*c).len() {
                        let snip = s.clone();
                        let m = Marker::Deletion;
                        let p = snip.pos as i32 + read_offset;
                        let f = snip.flags;
                        let n = snip.name;
                        let rs = snip.pos;
                        let re = snip.pos as i32 + snip.length as i32;
                        let b = String::from("");

                        let base = AlignmentNucleobase {
                            marker_type: m,
                            bases: b,
                            position: p as f32,
                            flags: f,
                            name: n,
                            read_start: rs as u32,
                            read_end: re as u32,
                        };

                        read_offset += 1;

                        if from as f32 <= base.position && base.position <= to as f32 {
                            bases.push(base);
                        }
                    }

                    soft_clip_begin = false;

                }
                rust_htslib::bam::record::Cigar::RefSkip(c) => {
                    for _i in 0..rust_htslib::bam::record::Cigar::RefSkip(*c).len() {
                        //offset += 1;
                    }

                    soft_clip_begin = false;

                }
                rust_htslib::bam::record::Cigar::SoftClip(c) => {
                    if soft_clip_begin {


                        for i in 0..rust_htslib::bam::record::Cigar::SoftClip(*c).len() {
                            let snip = s.clone();
                            let b = char_vec[cigar_offset as usize];

                            let dif = rust_htslib::bam::record::Cigar::SoftClip(*c).len() - i;

                            if snip.pos - dif as i32 >= from as i32 && (snip.pos - dif as i32) < to as i32 {
                                let ref_index = snip.pos - dif as i32 - from as i32;
                                let ref_base = &ref_bases[ref_index as usize];
                                let m: Marker;
                                if ref_base.get_marker_type() == b {
                                    m = Marker::Match; // Match with reference fasta
                                } else {
                                    match b { // Mismatch
                                        'A' => m = Marker::A,
                                        'T' => m = Marker::T,
                                        'C' => m = Marker::C,
                                        'N' => m = Marker::N,
                                        'G' => m = Marker::G,
                                        _ => m = Marker::Deletion,
                                    }
                                }




                                let p = snip.pos as i32 - dif as i32;
                                let f = snip.flags;
                                let n = snip.name;

                                let rs: i32;
                                let re: i32;

                                if snip.paired {
                                    if snip.pos < snip.mate_pos {
                                        re = snip.mate_pos + 100;
                                        rs = snip.pos - rust_htslib::bam::record::Cigar::SoftClip(*c).len() as i32;
                                    } else {
                                        rs = snip.mate_pos;
                                        re = snip.pos as i32 + snip.length as i32;
                                    }
                                } else {
                                    rs = snip.pos - rust_htslib::bam::record::Cigar::SoftClip(*c).len() as i32;
                                    re = snip.pos as i32 + snip.length as i32;
                                }


                                let base = AlignmentNucleobase {
                                    marker_type: m,
                                    bases: b.to_string(),
                                    position: p as f32,
                                    flags: f,
                                    name: n,
                                    read_start: rs as u32,
                                    read_end: re as u32,
                                };


                                bases.push(base);
                            }
                            cigar_offset += 1;
                        }
                    } else {
                        for _i in 0..rust_htslib::bam::record::Cigar::SoftClip(*c).len() {
                            let snip = s.clone();
                            let b = char_vec[cigar_offset as usize];

                            if snip.pos + read_offset >= from as i32 && snip.pos + read_offset < to as i32 {
                                let ref_index = snip.pos + read_offset - from as i32;
                                let ref_base = &ref_bases[ref_index as usize];


                                let m: Marker;
                                if ref_base.get_marker_type() == b {
                                    m = Marker::Match; // Match with reference fasta
                                } else {
                                    match b { // Mismatch
                                        'A' => m = Marker::A,
                                        'T' => m = Marker::T,
                                        'C' => m = Marker::C,
                                        'N' => m = Marker::N,
                                        'G' => m = Marker::G,
                                        _ => m = Marker::Deletion,
                                    }
                                }


                                let p = snip.pos as i32 + read_offset;
                                let f = snip.flags;
                                let n = snip.name;

                                let rs: i32;
                                let re: i32;

                                if snip.paired {
                                    if snip.pos < snip.mate_pos {
                                        re = snip.mate_pos + 100;
                                        rs = snip.pos;
                                    } else {
                                        rs = snip.mate_pos;
                                        re = snip.pos as i32 + snip.length as i32;
                                    }
                                } else {
                                    rs = snip.pos;
                                    re = snip.pos as i32 + snip.length as i32;
                                }


                                let base = AlignmentNucleobase {
                                    marker_type: m,
                                    bases: b.to_string(),
                                    position: p as f32,
                                    flags: f,
                                    name: n,
                                    read_start: rs as u32,
                                    read_end: re as u32,
                                };


                                bases.push(base);
                            }
                            cigar_offset += 1;
                            read_offset += 1;
                        }
                    }

                    soft_clip_begin = false;

                }
                rust_htslib::bam::record::Cigar::HardClip(c) => {
                    for _i in 0..rust_htslib::bam::record::Cigar::HardClip(*c).len() {
                        cigar_offset += 1;
                    }

                    soft_clip_begin = false;

                }
                rust_htslib::bam::record::Cigar::Pad(c) => {
                    for _i in 0..rust_htslib::bam::record::Cigar::Pad(*c).len() {
                        //offset += 1;
                    }

                    soft_clip_begin = false;

                }
                rust_htslib::bam::record::Cigar::Equal(c) => {
                    for _i in 0..rust_htslib::bam::record::Cigar::Equal(*c).len() {
                        //offset += 1;
                    }

                    soft_clip_begin = false;

                }
                rust_htslib::bam::record::Cigar::Diff(c) => {
                    for _i in 0..rust_htslib::bam::record::Cigar::Diff(*c).len() {
                        //offset += 1;
                    }

                    soft_clip_begin = false;

                }
            }
        }
        //println!("Cigar Offset:{}, Read Offset {}", cigar_offset, read_offset);
    }
    bases
}

pub fn get_reads(path: &Path, fasta_path: &Path, chrom: String, from: u32, to: u32) -> Vec<AlignmentNucleobase> {
    let alignments = read_indexed_bam(path,chrom.clone(), from, to);
    let bases = make_nucleobases(fasta_path, chrom, alignments, from, to);

    bases
}