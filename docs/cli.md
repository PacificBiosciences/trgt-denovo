## TRGT-denovo command-line options

### Command: trio and duo

Both the trio and duo commands are used to process tandem repeats, with trio designed for analyzing data from a mother, father, and child, and duo intended for cases with only two individuals. Below are the common and specific options for each command:

Basic options (Common to both trio and duo):

- `-r, --reference <FASTA>` Path to the FASTA file containing the reference genome. Use the same reference genome as for read alignment. Its index must be at `<FASTA>.fai`; create it with `samtools faidx <FASTA>`. The FASTA filename need not have an extension.
- `-b, --bed <BED>` Path to the BED file with reference coordinates of tandem repeats
- `-o, --out <TSV>` Output TSV path. Required for whole-catalog runs; optional with `--trid`, where it overrides stdout.
- `--alignment-scores <PATH>` Optionally export individual read-to-candidate alignment scores to a separate TSV file. Available in both modes and independent of `--verbose` and `--quiet`; see [alignment-score output](output.md#individual-read-alignment-scores). Use a different path from `--out`.
- `--trid <TRID>` Analyze only the first matching repeat ID in the BED file. Results go to stdout by default, or to the file specified by `--out`; default = None.
- `-@ <THREADS>` Number of threads, the number of sites that are processed in parallel, default = 1
- `-h, --help` Print help
- `-V, --version` Print version; available at the root command and under either subcommand.
- `-v/--verbose` Increase verbosity (repeat flag for more detail); `--quiet` silences diagnostic logging, not explicitly requested output files.

Input file paths must name regular files. The parent directories of sample prefixes and any explicit output path must already exist, and the output path cannot be a directory. BAM and VCF files must also have valid indexes; format and index errors are reported when the files are opened.

Options specific to trio:

- For each sample you must supply either a prefix or explicit files (mixing across samples is allowed):
  - `-m, --mother <PREFIX>` or `--mother-vcf <VCF>` **and** `--mother-bam <BAM>`
  - `-f, --father <PREFIX>` or `--father-vcf <VCF>` **and** `--father-bam <BAM>`
  - `-c, --child <PREFIX>` or `--child-vcf <VCF>` **and** `--child-bam <BAM>`

Options specific to duo:

- For each sample you must supply either a prefix or explicit files:
  - `-1, --sample-a <PREFIX>` or `--sample-a-vcf <VCF>` **and** `--sample-a-bam <BAM>`
  - `-2, --sample-b <PREFIX>` or `--sample-b-vcf <VCF>` **and** `--sample-b-bam <BAM>`

Advanced:

- `--flank-len <FLANK_LEN>` Amount of additional flanking sequence that should be used during alignment, default = 50
- `--no-clip-aln` Score alignments without stripping the flanks
- `--p-quantile <QUANTILE>` Quantile of alignment scores to determine the threshold, default is strict and takes only the top scoring alignment, default = 1.0
- `--aln-scoring <SCORING>` Two-piece affine alignment penalties in the order `mismatch,gap_opening1,gap_extension1,gap_opening2,gap_extension2`, default = `8,4,2,24,1`. Supply exactly five integers: mismatch and gap extensions must be positive, while gap openings may be zero. Empty, non-integer, extra, and out-of-range integer fields are rejected. The supplied penalties apply to all alignments in either command, including multithreaded runs. See the [WFA2-lib repository](https://github.com/smarco/WFA2-lib/) for more details on parametrization.
- `--partition-by-aln` Within-sample partitioning using alignment rather than the TRGT BAMleet allele length field.
- `--quick <QUICK>` Only test loci that differ by a certain fraction in allele length. Format: `<field>,<fraction>` (e.g. AL,0.1 or AL). If no fraction is specified (or fraction is 0), it checks for exact matches. If a fraction is specified, it checks if the relative difference is within the given tolerance

Verbose output:

You can increase the verbosity of TRGT-denovo's output by adding a verbose flag before the command (e.g., `trgt-denovo -vv`):

- `-v` Provides verbose output
- `-vv` Provides even more detailed verbose output
