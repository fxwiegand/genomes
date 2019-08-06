use bio::io::fasta;
use std::path::Path;

pub fn read_fasta(path: &Path, chrom: String, start: u64, stop: u64) -> Vec<Nucleobase>    {
    let mut reader = fasta::IndexedReader::from_file(&path).unwrap();
    let index = fasta::Index::with_fasta_file(&path).unwrap();
    let _sequences = index.sequences();

    let mut seq:Vec<u8> = Vec::new();

    //println!("Reading genome number {}.", &seq_name.name);

    reader.fetch(&chrom, start, stop).unwrap();
    reader.read(& mut seq);

    let mut fasta = Vec::new();
    let mut ind = start;
    for a in seq {
        let b = Nucleobase {
            start_position: ind as f64 - 0.5,
            end_position: ind as f64 + 0.5,
            marker_type: a as char,
            row: 0
        };
        fasta.push(b);
        ind += 1;
    }

    fasta

}

pub fn read_fasta_at_pos(path: &Path, chrom: String, pos: u64) -> Nucleobase {
    let mut reader = fasta::IndexedReader::from_file(&path).unwrap();
    let index = fasta::Index::with_fasta_file(&path).unwrap();
    let _sequences = index.sequences();

    let mut seq:Vec<u8> = Vec::new();

    reader.fetch(&chrom, pos, pos + 1).unwrap();
    reader.read(& mut seq);


    let a = seq[0];
    let base = Nucleobase {
        start_position: pos as f64 - 0.5,
        end_position: pos as f64 + 0.5,
        marker_type: a as char,
        row: 0
    };

    base
}

#[derive(Serialize, Clone)]
pub struct Nucleobase {
    start_position: f64,
    end_position: f64,
    marker_type: char,
    row: u8,
}

impl Nucleobase {
    pub fn get_marker_type(&self) -> char { self.marker_type }
}