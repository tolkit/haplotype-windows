pub mod var_types {
    pub struct VarTypes {
        pub snps: i32,
        pub insertions: i32,
        pub deletions: i32,
        pub snp_density: i32,
        pub transitions: i32,
        pub transversions: i32,
        pub mean_indel_length: f32,
    }

    fn filter_and_average(numbers: Vec<isize>) -> f32 {
        let filt_nums: Vec<&isize> = numbers
            .iter()
            .filter(|&i| *i > 0isize)
            .collect::<Vec<&isize>>();

        filt_nums.clone().into_iter().sum::<isize>() as f32 / filt_nums.clone().len() as f32
    }

    pub fn count_var_types(alleles: &mut Vec<(Vec<u8>, Vec<u8>)>) -> VarTypes {
        // variant type
        let mut snps = 0;
        let mut insertions = 0;
        let mut deletions = 0;
        // transitions / transversions
        let mut transitions = 0;
        let mut transversions = 0;
        let mut indel_lengths: Vec<isize> = Vec::new();

        // iterate over the tuples in 'alleles' (ref, alt)
        for (a, b) in alleles.clone() {
            // collect a vector of absolute indel lengths
            // TODO: discuss whether this actually makes sense
            let ia_len: isize = a.len() as isize;
            let ib_len: isize = b.len() as isize;

            indel_lengths.push(ia_len.wrapping_sub(ib_len).abs());

            // number of snps, insertions, deletions
            if a.len() != 0 && b.len() != 0 {
                if a.len() == b.len() {
                    snps += 1;
                } else if b.len() > a.len() {
                    deletions += 1;
                } else if a.len() > b.len() {
                    insertions += 1;
                }
            }
            // number of transitions and transversions
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
        let mean_indel_length = filter_and_average(indel_lengths);

        VarTypes {
            snps,
            insertions,
            deletions,
            snp_density,
            transitions,
            transversions,
            mean_indel_length,
        }
    }
}
