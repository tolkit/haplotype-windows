pub mod windows {

    // this module writes to file all the calcualted stats
    // TODO: is there a more elegant way to do this?

    use std::io::LineWriter;
    use std::io::Write;

    pub fn window_calcs<T: std::io::Write>(current_contig: &str, current_record: &rust_htslib::bcf::Record, current_rid: &mut u32, current_window: &mut i32, window_size: usize, file2: &mut LineWriter<T>, alleles: &mut Vec<(Vec<u8>, Vec<u8>)>) {

            if current_record.rid().unwrap() == *current_rid {

                // this code chunk checks whether there is a gap between variants which is larger than the window size
                // if so, fill these gaps with zero records.
                if (current_record.pos() as i32) - *current_window > window_size as i32 {
                    let x = (current_record.pos() as i32) - *current_window;
                    let iterations = (x as f64 / window_size as f64).ceil() as i32;

                    if iterations > 0 {
                        for _it in 0..iterations {
                            writeln!(file2, "{},{},{},{},{},{},{},{}", current_contig, current_window, 0, 0, 0, 0, 0, 0).unwrap_or_else(|_| println!("[-]\tError in writing to file."));
                            *current_window += window_size as i32;
                        }
                    }

                }
                // collect records up until the current window
                if (current_record.pos() as i32) <= *current_window {

                    // collect the window records
                    let reference = &current_record.alleles()[0].to_vec();
                    let alternate = &current_record.alleles()[1].to_vec();
                    
                    // now, alleles is a vector of u8's
                    alleles.push((reference.to_vec(), alternate.to_vec()));
    
                // at the record which > current window, process the records
                } else if (current_record.pos() as i32) > *current_window { 

                    // variant type
                    let mut snps = 0;
                    let mut insertions = 0;
                    let mut deletions = 0;
                    // transitions / transversions
                    let mut transitions = 0;
                    let mut transversions = 0;
    
                    for (a, b) in alleles.clone() {
                        if a.len() != 0 && b.len() != 0 {
                            if a.len() == b.len() {
                                snps += 1;
                            } else if b.len() > a.len() {
                                deletions += 1;
                            } else if a.len() > b.len() {
                                insertions += 1;
                            }
                        }
                        // transitions
                        // A G
                        if a == [65] && b == [71] || a == [71] && b == [65] {
                            transitions += 1;
                        } else 
                        // C T
                        if a == [67] && b == [84] || a == [84] && b == [67] {
                            transitions += 1;
                        } else 
                        // transversions
                        // A C
                        if a == [65] && b == [67] || a == [67] && b == [65] {
                            transversions += 1;
                        } else
                        // G T
                        if a == [71] && b == [84] || a == [84] && b == [71] {
                            transversions += 1;
                        } else 
                        // A T
                        if a == [65] && b == [84] || a == [84] && b == [65] {
                            transversions += 1;
                        } else 
                        // C G
                        if a == [67] && b == [71] || a == [71] && b == [67] {
                            transversions += 1;
                        }
                    }
                    // total number of SNPs
                    let snp_density = snps + deletions + insertions;
                    
                    if alleles.is_empty() {
                        println!("[-]\tIf this message appears, there are no alleles in the current position of the genome. This is clearly a bug!");
                    } else {
                        writeln!(file2, "{},{},{},{},{},{},{},{}", current_contig, current_window, snps, insertions, deletions, snp_density, transitions, transversions).unwrap_or_else(|_| println!("-]\tError in writing to file."));
                    }
    
                    // clear the collections
                    alleles.clear();
                    // add the alleles of the skipped row... [is this right?]
                    let reference = &current_record.alleles()[0].to_vec();
                    let alternate = &current_record.alleles()[1].to_vec();
                    alleles.push((reference.to_vec(), alternate.to_vec()));
                    // else, if the position in the genome is greater than the current window, then increase the size of the window
                    *current_window += window_size as i32;
                }
            } else {
                // go to the next record and reset the window
                *current_rid += 1;
                *current_window = window_size as i32;
            }
    }
}