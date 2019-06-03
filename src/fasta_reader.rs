use bio::io::fasta;
use std::path::Path;

pub fn read_fasta(path: &Path, chrom: u8, start: u64, stop: u64) -> Vec<Nucleobase>    {
    let mut reader = fasta::IndexedReader::from_file(&path).unwrap();
    let index = fasta::Index::with_fasta_file(&path).unwrap();
    let sequences = index.sequences();
    let seq_name = &sequences[chrom as usize];

    let mut seq:Vec<u8> = Vec::new();

    println!("Reading genome number {}.", &seq_name.name);

    reader.fetch(&seq_name.name, start, stop).unwrap();
    reader.read(& mut seq);

    let mut fasta = Vec::new();
    let mut ind = start;
    for a in seq {
        let b = Nucleobase {
            position: ind,
            base: a as char,
            row: 0
        };
        fasta.push(b);
        ind += 1;
    }

    fasta

}

#[derive(Serialize, Clone)]
pub struct Nucleobase {
    position: u64,
    base: char,
    row: u8,
}