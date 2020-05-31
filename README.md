![build & test](https://github.com/koesterlab/genomes/workflows/build%20&%20test/badge.svg)

# Genomes

A web based tool for genome visualization that provides these features:

* A full client-server application for visualizing genome data (`cargo run server`)
* Creating static vega plots (`cargo run static`)
* Creating a html report for VCF files (`cargo run report`)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

Don't forget to index your fasta-, bam- and vcf-files with samtools

```
samtools index path/to/myfasta.fa 
```

### Running

Start the server with:

```
cargo run server data/mybam.bam data/myfasta.fa data/myvcf.vcf.gz
```

For a static view that you can paste into the [Vega Online Editor](https://vega.github.io/editor/) or render with the [Vega Command Line Utilities](https://vega.github.io/vega/usage/#cli) start with:

```
cargo run static data/mybam.bam data/myfasta.fa data/myvcf.vcf.gz chromosom from to > visualization.json
```

To create a html report run the following:

```
cargo run report -r data/mybam.bam data/myfasta.fa data/myvcf.vcf.gz chromosom > report.html
```
or without the `-r` flag to start a server that deploys the html as a website on your local machine

## Built With

* [Rocket](https://rocket.rs) - A web framework for Rust
* [Vega/Vega-Lite](https://vega.github.io) - A visualization grammar
* [Rust-Htslib](https://github.com/rust-bio/rust-htslib) - HTSlib bindings and a high level Rust API for reading and writing BAM files
* [Rust-Bio](https://github.com/rust-bio/rust-bio) - algorithms and data structures that are useful for bioinformatics
* [Tera](https://tera.netlify.app) - A powerful, easy to use template engine for Rust

## Authors

* **Felix Wiegand** - (https://github.com/fxwiegand)

See also the list of contributors who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Created for the bachelor thesis of **Felix Wiegand** as TU Dortmund University


