pub mod density {

    use std::collections::HashMap;
    use std::io::LineWriter;
    use itertools::Itertools;
    use std::io::Write;

    pub fn snp_density<T: std::io::Write>(current_contig: &str, current_record: &rust_htslib::bcf::Record, current_rid: &mut u32, current_window: &mut i32, window_size: usize, map: &mut HashMap<i32, i32>, file: &mut LineWriter<T>) {

        if current_record.rid().unwrap() == *current_rid {
            // if the current window is less than the position in the genome, increment
            if (current_record.pos() as i32) < *current_window {
                let count = map.entry(*current_window).or_insert(0);
                *count += 1;
            } else { // else, increase the size of the window (is this right?)
                *current_window += window_size as i32;
                // print
                for key in map.keys().sorted() {
                    writeln!(file, "{},{},{}", current_contig, key, map[key]).unwrap();
                }
                // clear and increment
                map.clear();
            }
        } else {
            *current_rid += 1;
            *current_window = 0;
        }
    }
}