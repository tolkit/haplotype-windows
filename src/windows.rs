pub mod windows {

    // this module writes to file all the calcualted stats
    // TODO: is there a more elegant way to do this?

    use std::io::LineWriter;
    use std::io::Write;
    //use rust_htslib::bcf::header::Id;

    pub fn window_calcs<T: std::io::Write>(current_contig: &str, current_record: &rust_htslib::bcf::Record, current_rid: &mut u32, current_window: &mut i32, window_size: usize, file2: &mut LineWriter<T>, alleles: &mut Vec<(i32, i32)>) {

        // TODO: deepvariant has this 'RefCall' on the filter, do we ignore these?
        // if we need to deal with this, wrap in: if current_record.has_filter(Id(0)) {...}
            if current_record.rid().unwrap() == *current_rid {

                // watch out, this is some hacky bullshit.
                if (current_record.pos() as i32) - *current_window > window_size as i32 {
                    let x = (current_record.pos() as i32) - *current_window;
                    let iterations = (x as f64 / window_size as f64).ceil() as i32;

                    if iterations > 0 {
                        for _it in 0..iterations {
                            writeln!(file2, "{},{},{},{},{},{}", current_contig, current_window, 0, 0, 0, 0).unwrap_or_else(|_| println!("[-]\tError in writing to file."));
                            *current_window += window_size as i32;
                        }
                    }

                }

                if (current_record.pos() as i32) <= *current_window {
    
                    // SNP vs INDELS
                    // collect the window records here *
                    let reference = &current_record.alleles()[0].to_vec();
                    let alternate = &current_record.alleles()[1].to_vec();
    
                    // lengths indicate indels...
                    alleles.push((reference.to_vec().len() as i32, alternate.to_vec().len() as i32));
    
                } else if (current_record.pos() as i32) > *current_window { 

                    // write snp type to file
                    let mut snps = 0;
                    let mut insertions = 0;
                    let mut deletions = 0;
    
                    for (a, b) in alleles.clone() {
                        if a != 0 && b != 0 {
                            if a == b {
                                snps += 1;
                            } else if b > a {
                                deletions += 1;
                            } else if a > b {
                                insertions += 1;
                            }
                        }
                    }
                    // total number of SNPs
                    let snp_density = snps + deletions + insertions;
                    
                    if alleles.is_empty() {
                        println!("[-]\tIf this message appears, there are no alleles in the current position of the genome. This is clearly a bug!");
                    } else {
                        writeln!(file2, "{},{},{},{},{},{}", current_contig, current_window, snps, insertions, deletions, snp_density).unwrap_or_else(|_| println!("-]\tError in writing to file."));
                    }
    
                    // clear and increment
                    alleles.clear();
                    // add the alleles of the skipped row... is this right??
                    let reference = &current_record.alleles()[0].to_vec();
                    let alternate = &current_record.alleles()[1].to_vec();
                    alleles.push((reference.to_vec().len() as i32, alternate.to_vec().len() as i32));
                    // else, if the position in the genome is greater than the current window, then increase the size of the window
                    *current_window += window_size as i32;
                }
            } else {
                // go to the next record!
                *current_rid += 1;
                *current_window = window_size as i32;
            }
    }
}