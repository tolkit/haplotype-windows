// standard library
use std::fs::{File, create_dir_all};
use std::io::LineWriter;
use std::str;
use std::io::prelude::*;

// non standard
use rust_htslib::bcf::{Reader, Read};
use clap::{App, Arg, value_t};

// my modules
use haplotype_windows::windows::windows;

// TODOs
// add number of substitution types?
// e.g. A <> C (transversion), C <> T (transition) etc.


fn main() {
    // command line options
    let matches = App::new("Haplotype windows")
        .version(clap::crate_version!())
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("Haplotype stats in windows on a VCF file.")
        .arg(Arg::with_name("vcf")
                 .short("f")
                 .long("vcf")
                 .takes_value(true)
                 .required(true)
                 .help("The input vcf (or bcf) file. Gzipped or not."))
        .arg(Arg::with_name("window_size")
                 .short("w")
                 .long("window_size")
                 .help("Integer size of window for statistics to be computed over.")
                 .takes_value(true)
                 .default_value("10000"))
        .arg(Arg::with_name("output")
                 .short("o")
                 .long("output")
                 .help("Output filename for the CSVs (without extension).")
                 .takes_value(true)
                 .required(true))
        .get_matches();

    // parse command line options
    let input_vcf = matches.value_of("vcf").unwrap();
    let output = matches.value_of("output").unwrap();
    let window_size = value_t!(matches.value_of("window_size"), usize).unwrap_or_else(|e| e.exit());

    // create directory for output
    if let Err(e) = create_dir_all("./hw_out/") {
        println!("[-]\tCreate directory error: {}", e.to_string());   
    }

    // initiate the output CSV for SNP type
    let output_file_2 = format!("./hw_out/{}{}", output, "_snp_type.csv");
    let window_file_2 = File::create(&output_file_2).unwrap();
    let mut window_file_2 = LineWriter::new(window_file_2);
    // add headers
    writeln!(window_file_2, "ID,window,no_snps,no_insertions,no_deletions,snp_density").unwrap();

    // the current window, incremented
    let mut current_window: i32 = window_size as i32;
    // i32 representation of the chromosomes to iterate over
    let mut curr_rid = 0;
    // where we collect the alleles repeatedly in the windows::window_calcs function below
    let mut window_alleles = Vec::new();

    // read in the VCF/BCF
    let mut vcf = Reader::from_path(input_vcf).expect("Error opening file.");
    
    // iterate over the VCF records
    for record_result in vcf.records() {
        let record = record_result.expect("Fail to read record");

        // u8 -> str, so a human can (easily) read the chromosome names
        let contig = record.header().rid2name(record.rid().unwrap()).unwrap();
        let s = match str::from_utf8(contig) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };
        windows::window_calcs(s, &record, &mut curr_rid, &mut current_window, window_size, &mut window_file_2, &mut window_alleles);
    }
}