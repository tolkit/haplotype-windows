use std::fs::{File, create_dir_all};
use std::io::LineWriter;
use std::str;
use std::collections::HashMap;
use std::io::prelude::*;

use itertools::Itertools;
use rust_htslib::bcf::{Reader, Read};
use clap::{App, Arg, value_t};

// TODOs
// add number of substitution types?
// e.g. A <> C (transversion), C <> T (transition) etc.
// types of variants e.g. SNPs/indels...
// can I convert the window generation to a function?


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
                 .help("The input vcf (or bcf) file."))
        .arg(Arg::with_name("window_size")
                 .short("w")
                 .long("window_size")
                 .help("Integer size of window for statistics to be computed over.")
                 .takes_value(true)
                 .default_value("10000"))
        .arg(Arg::with_name("output")
                 .short("o")
                 .long("output")
                 .help("Output filename for the CSV (without extension).")
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
    // initiate the output CSV for windows
    let output_file_1 = format!("./hw_out/{}{}", output, "_windows.csv");
    let window_file = File::create(&output_file_1).unwrap();
    let mut window_file = LineWriter::new(window_file);
    // add headers
    writeln!(window_file, "ID,window,SNP_density").unwrap();

    let mut window_snp_map: HashMap<i32, i32> = HashMap::new(); // the window and number of SNPs
    let mut current_window: i32 = window_size as i32;
    let mut curr_rid = 0;

    let mut bcf = Reader::from_path(input_vcf).expect("Error opening file.");
    
    for record_result in bcf.records() {
        let record = record_result.expect("Fail to read record");
        let contig = record.header().rid2name(record.rid().unwrap()).unwrap(); // ... made it hard to get the contig name back out!
        let s = match str::from_utf8(contig) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };

        if record.rid().unwrap() == curr_rid {
            // if the current window is less than the position in the genome, increment
            if (record.pos() as i32) < current_window {
                let count = window_snp_map.entry(current_window).or_insert(0);
                *count += 1;
            } else { // else, increase the size of the window (is this right?)
                current_window += window_size as i32;
                // print
                for key in window_snp_map.keys().sorted() {
                    writeln!(window_file, "{},{},{}", s, key, window_snp_map[key]).unwrap();
                }
                // clear and increment
                window_snp_map.clear();
            }
        } else {
            curr_rid += 1;
            current_window = 0;
        }
    }
}