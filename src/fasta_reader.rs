use bio::io::fasta;
use std::path::Path;

pub fn read_fasta(path: &Path, chrom: u8, start: u64, stop: u64) -> String    {
    let mut reader = fasta::IndexedReader::from_file(&path).unwrap();
    let index = fasta::Index::with_fasta_file(&path).unwrap();
    let sequences = index.sequences();
    let seq_name = &sequences[chrom as usize];

    let mut seq:Vec<u8> = Vec::new();

    println!("Reading genome number {}.", &seq_name.name);

    reader.fetch(&seq_name.name, start, stop).unwrap();
    reader.read(& mut seq);

    let mut fasta = String::from("");
    for a in seq {
        fasta.push(a as char);
    }

    fasta
}


