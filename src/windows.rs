pub mod windows {

    // this module writes to file all the calcualted stats
    // TODO: is there a more elegant way to do this?

    // std modules
    use std::io::LineWriter;
    use std::io::Write;
    // my modules
    use crate::var_types::var_types;

    pub fn window_calcs<T: std::io::Write>(
        current_contig: &str,
        current_record: &rust_htslib::bcf::Record,
        current_rid: &mut u32,
        current_window: &mut i32,
        window_size: usize,
        file2: &mut LineWriter<T>,
        alleles: &mut Vec<(Vec<u8>, Vec<u8>)>,
    ) {
        if current_record.rid().unwrap() == *current_rid {
            // this code chunk checks whether there is a gap between variants which is larger than the window size
            // if so, fill these gaps with zero records.
            if (current_record.pos() as i32) - *current_window > window_size as i32 {
                let x = (current_record.pos() as i32) - *current_window;
                let iterations = (x as f64 / window_size as f64).ceil() as i32;

                if iterations > 0 {
                    for _it in 0..iterations {
                        writeln!(
                            file2,
                            "{},{},{},{},{},{},{},{},{}",
                            current_contig, current_window, 0, 0, 0, 0, 0, 0, 0
                        )
                        .unwrap_or_else(|_| println!("[-]\tError in writing to file."));
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
                // count the different kinds of variant.
                let var_counts = var_types::count_var_types(alleles);

                if alleles.is_empty() {
                    // I think these are invariant sites?
                    println!("[-]\tInvariant site found.");
                } else {
                    writeln!(
                        file2,
                        "{},{},{},{},{},{},{},{},{}",
                        current_contig,
                        current_window,
                        var_counts.snps,
                        var_counts.insertions,
                        var_counts.deletions,
                        var_counts.snp_density,
                        var_counts.transitions,
                        var_counts.transversions,
                        var_counts.mean_indel_length
                    )
                    .unwrap_or_else(|_| println!("[-]\tError in writing to file."));
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
