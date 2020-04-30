# Genomes

A web based tool for genome visualization 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

Don't forget to index your fasta-, bam- and vcf-files with samtools

```
samtools index path/to/myfasta.fa 
```

### Running

Start with:

```
cargo run server data/mybam.bam data/myfasta.fa data/myvcf.vcf.gz
```

For a static view that you can paste into the [Vega Online Editor](https://vega.github.io/editor/) or render with the [Vega Command Line Utilities](https://vega.github.io/vega/usage/#cli) start with:

```
cargo run static data/mybam.bam data/myfasta.fa data/myvcf.vcf.gz chromosom from to > visualization.json
```

## Built With

* [Rocket](https://rocket.rs) - A web framework for Rust
* [Vega/Vega-Lite](https://vega.github.io) - A visualization grammar
* [Rust-Htslib](https://github.com/rust-bio/rust-htslib) - HTSlib bindings and a high level Rust API for reading and writing BAM files
* [Rust-Bio](https://github.com/rust-bio/rust-bio) - algorithms and data structures that are useful for bioinformatics

## Authors

* **Felix Wiegand** - (https://github.com/fxwiegand)

See also the list of contributors who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Created for the bachelor thesis of **Felix Wiegand** as TU Dortmund University


