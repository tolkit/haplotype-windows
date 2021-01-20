# Comparing DToL alternate haplotypes

## Scope

Working on a VCF (or BCF, gzipped or not), calculate SNP density in windows. Very much in progress.

## TODOs

- Perhaps count types of variant (indels/SNPs)
- Count transversions/transitions
- Anything else?

## Current usage

You have to complile yourself. <a href="https://www.rust-lang.org/tools/install">Download rust</a>, clone this repo, and then run:

`cargo build --release`

Compiling may take a couple of minutes. This will then make the compiled binary in the `target/release` directory.

Run `./target/release/fasta_windows` --help to display the help message in the terminal.

```
Haplotype windows 0.1.0
Max Brown <mb39@sanger.ac.uk>
Haplotype stats in windows on a VCF file.

USAGE:
    haplotype_windows [OPTIONS] --output <output> --vcf <vcf>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -o, --output <output>              Output filename for the CSV (without extension).
    -f, --vcf <vcf>                    The input vcf (or bcf) file.
    -w, --window_size <window_size>    Integer size of window for statistics to be computed over. [default: 10000]
```

Outputs a CSV:

```
ID,window,SNP_density
SUPER_1,10000,42
SUPER_1,20000,32
SUPER_1,30000,53
SUPER_1,40000,195
SUPER_1,50000,245
SUPER_1,60000,158
SUPER_1,70000,74
SUPER_1,80000,225
SUPER_1,90000,202
SUPER_1,100000,325
```
