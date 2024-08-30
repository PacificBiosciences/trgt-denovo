## TRGT-denovo command-line options

### Command: trio and duo

Both the trio and duo commands are used to process tandem repeats, with trio designed for analyzing data from a mother, father, and child, and duo intended for cases with only two individuals. Below are the common and specific options for each command:

Basic Options (Common to both trio and duo):
- `-r, --reference <FASTA>` Path to the FASTA file containing the reference genome. The same reference genome as was used for read alignment.
- `-b, --bed <BED>` Path to the BED file with reference coordinates of tandem repeats
- `-o, --out <TSV>` Output tsv path
- `--trid <TRID>` Optionally a tandem repeat ID may be supplied, if so, only this site will be tested, default = None
- `-@ <THREADS>` Number of threads, the number of sites that are processed in parallel, default = 1
- `-h, --help` Print help
- `-V, --version` Print version

Options specific to trio:
- `-m, --mother <PREFIX>` Common (path) prefix of spanning reads BAM file and variant call VCF file of mother
- `-f, --father <PREFIX>` Common (path) prefix of spanning reads BAM file and variant call VCF file of father
- `-c, --child <PREFIX>` Common (path) prefix of spanning reads BAM file and variant call VCF file of child

Options specific to duo:
- `-a, --sample-a <PREFIX>`  Common (path) prefix of spanning reads BAM file and variant call VCF file of the first sample
- `-b, --sample-b <PREFIX>`  Common (path) prefix of spanning reads BAM file and variant call VCF file of the second sample

Advanced:
- `--flank-len <FLANK_LEN>` Amount of additional flanking sequence that should be used during alignment, default = 50
- `--no-clip-aln` Score alignments without stripping the flanks
- `--p-quantile <QUANTILE>` Quantile of alignment scores to determine the threshold, default is strict and takes only the top scoring alignment, default = 1.0
- `--aln-scoring` Scoring function for 2-piece gap affine alignment (non-negative values): mismatch,gap_opening1,gap_extension1,gap_opening2,gap_extension2, default = "8,4,2,24,1", see [here](https://github.com/smarco/WFA2-lib/) for more details on parametrization
- `--partition-by-aln` Within-sample partitioning using alignment rather than the TRGT BAMleet allele length field.
- `--quick <QUICK>` Only test loci that differ by a certain fraction in allele length. Format: <field>,<fraction> (e.g. AL,0.1 or AL). If no fraction is specified (or fraction is 0), it checks for exact matches. If a fraction is specified, it checks if the relative difference is within the given tolerance.

Verbose Output:

You can increase the verbosity of TRGT-denovo's output by adding a verbose flag before the command:
- `-v` Provides verbose output.
- `-vv` Provides more detailed verbose output.