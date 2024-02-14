## TRGT-denovo command-line options

### Command: trio
Basic:
- `-r, --reference <FASTA>` Path to the FASTA file containing the reference genome. The same reference genome as was used for read alignment.
- `-b, --bed <BED>` Path to the BED file with reference coordinates of tandem repeats
- `-m, --mother <PREFIX>` Common (path) prefix of spanning reads BAM file and variant call VCF file of mother
- `-f, --father <PREFIX>` Common (path) prefix of spanning reads BAM file and variant call VCF file of father
- `-c, --child <PREFIX>` Common (path) prefix of spanning reads BAM file and variant call VCF file of child
- `-o, --out <CSV>` Output csv path
- `--trid <TRID>` Optionally a tandem repeat ID may be supplied, if so, only this site will be tested, default = None
- `-@ <THREADS>` Number of threads, the number of sites that are processed in parallel, default = 1
- `-h, --help` Print help
- `-V, --version` Print version

Advanced:
- `--flank-len <FLANK_LEN>` Number of flanking nucleotides added to a target region during realignment, default = 50
- `--no-clip-aln` Score alignments without stripping the flanks
- `--parental-quantile <QUANTILE>` Quantile of alignment scores to determine the parental threshold, default is strict and takes only the top scoring alignment, default = 1.0
- `--aln-scoring` Scoring function for 2-piece gap affine alignment (non-negative values): mismatch,gap_opening1,gap_extension1,gap_opening2,gap_extension2, default = "8,4,2,24,1", see [here](https://github.com/smarco/WFA2-lib/) for more details on parametrization 