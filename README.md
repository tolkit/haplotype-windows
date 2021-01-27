# Comparing DToL alternate haplotypes

## Scope

Working on a VCF (or BCF, gzipped or not), calculate SNP density in windows. Very much in progress. May be buggy.

I had a really (really) hard time thinking about how to split a VCF into windows (actually non-overlapping chunks), while passing the file only once. My approach definitely needs some quality control, but seems to match up to bedtools (albeit limited testing).

## Current usage

You have to complile yourself. <a href="https://www.rust-lang.org/tools/install">Download rust</a>, clone this repo, and then run:

`cargo build --release`

Compiling may take a couple of minutes. This will then make the compiled binary in the `target/release` directory.

Run `./target/release/fasta_windows` --help to display the help message in the terminal.

```
Haplotype windows 0.1.1
Max Brown <mb39@sanger.ac.uk>
Haplotype stats in windows on a VCF file.

USAGE:
    haplotype_windows [OPTIONS] --output <output> --vcf <vcf>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -o, --output <output>              Output filename for the CSVs (without extension).
    -p, --pass_only <pass_only>        Should calculations be output only on records with as PASS filter? Boolean, input
                                       true or false. [default: true]
    -f, --vcf <vcf>                    The input vcf (or bcf) file. Gzipped or not.
    -w, --window_size <window_size>    Integer size of window for statistics to be computed over. [default: 10000]
```

Outputs a CSV:

Currently spits out the numbers of SNPs, insertions, deletions, transitions, transversions, as well as the total number of variants per window.

e.g. `./target/release/haplotype-windows -f <vcf> -o test -p false`

```
ID,window,no_snps,no_insertions,no_deletions,snp_density,transitions,transversions
SUPER_1,10000,34,3,5,42,18,16
SUPER_1,20000,28,2,3,33,17,11
SUPER_1,30000,51,2,1,54,26,25
SUPER_1,40000,175,9,12,196,91,84
SUPER_1,50000,206,22,18,246,116,90
SUPER_1,60000,131,15,13,159,59,72
SUPER_1,70000,62,7,6,75,28,34
SUPER_1,80000,195,20,11,226,109,86
SUPER_1,90000,176,19,8,203,102,74
```
