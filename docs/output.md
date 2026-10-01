# Interpreting TRGT-denovo output

Missing values are denoted by: `.`, indicating that a given locus has a missing genotype in at least one sample, or a locus is skipped because of quick mode.

Trio fields:

- `chrom` contig of tandem repeat
- `start` start position of tandem repeat
- `end` end position of tandem repeat
- `motifs` motif definitions of tandem repeat
- `trid` ID of the tandem repeat, encoded as in the BED file
- `genotype` Genotype ID of the child for a specific allele, corresponds to the TRGT genotype ID
- `denovo_coverage` Number of child reads supporting a *de novo* allele compared to parental data
- `allele_coverage` Number of child reads mapped to this specific allele
- `allele_ratio` Ratio of *de novo* coverage to allele coverage
- `child_coverage` Total number of child reads at this site
- `child_ratio` Ratio of *de novo* coverage to total coverage at this site
- `mean_diff_father` Score difference between *de novo* and paternal reads; lower values indicate greater similarity
- `mean_diff_mother` Score difference between *de novo* and maternal reads; lower values indicate greater similarity
- `father_dropout_prob` Dropout rate for reads coming from the mother
- `mother_dropout_prob` Dropout rate for reads coming from the father.
- `allele_origin` Inferred origin of the allele based on alignment; possible values: `{F:{1,2,?}, M:{1,2,?}, ?}`. `F` and `M` denote father and mother respectively. The associated `{1, 2, ?}` values denote the first or second allele from either parent or `?` when this cannot be derived unambiguously. Lastly a `?` denotes an allele for which parental origin cannot be determined unambiguously
- `denovo_status` Indicates if the allele is *de novo*, only if `allele_origin` is defined; possible values: `{., X, Y:{+, -, =, ?}}`. This is `.` if there is a missing value, `X` if no *de novo* read is found and `Y` otherwise, if parental origin can be determined without ambiguity the allele sequences can be compared directly such that the *de novo* type can be established as `+` (expansion), `-` (contraction), `=` (substitution), `?` if `allele_origin` is not defined
- `per_allele_reads_father` Number of reads partitioned per allele in the father (allele1, allele2)
- `per_allele_reads_mother` Number of reads partitioned per allele in the mother (allele1, allele2)
- `per_allele_reads_child` Number of reads partitioned per allele in the child (allele1, allele2)
- `father_dropout` Coverage cut-off dropout detection using HP tags from phasing tools in father; possible values: Full dropout (`FD`), Haplotype dropout (`HD`), Not (`N`)
- `mother_dropout` Coverage cut-off dropout detection using HP tags from phasing tools in mother; possible values: Full dropout (`FD`), Haplotype dropout (`HD`), Not (`N`)
- `child_dropout` Coverage cut-off dropout detection using HP tags from phasing tools in child; possible values: Full dropout (`FD`), Haplotype dropout (`HD`), Not (`N`)
- `index` Index of this allele in the TRGT VCF
- `father_MC` TRGT VCF motif counts for this locus in the father
- `mother_MC` TRGT VCF motif counts for this locus in the mother
- `child_MC` TRGT VCF motif counts for this locus in the child
- `father_AL` TRGT VCF allele lengths for this locus in the father
- `mother_AL` TRGT VCF allele lengths for this locus in the mother
- `child_AL` TRGT VCF allele lengths for this locus in the child
- `father_overlap_coverage` Reciprocal of `denovo_coverage`, the number of reads in per allele in the father that overlap compared to the child data
- `mother_overlap_coverage` Reciprocal of `denovo_coverage`, the number of reads in per allele in the mother that overlap compared to the child data

Duo fields:

- `chrom` contig of tandem repeat
- `start` start position of tandem repeat
- `end` end position of tandem repeat
- `motifs` motif definitions of tandem repeat
- `trid` ID of the tandem repeat, encoded as in the BED file
- `genotype` Genotype ID of sample A for a specific allele, corresponds to the TRGT genotype ID
- `denovo_coverage` Number of sample A reads supporting a *de novo* allele compared to sample B
- `allele_coverage` Number of sample A reads mapped to this specific allele
- `allele_ratio` Ratio of *de novo* coverage to allele coverage
- `a_coverage` Total number of sample A reads at this site
- `a_ratio` Ratio of *de novo* coverage to total coverage at this site
- `mean_diff_b` Score difference between *de novo* and sample B reads; lower values indicate greater similarity
- `denovo_status` Attempts to say if the allele is *de novo*; possible values: `{., X, Y:{?}}`. This is `.` if there is a missing value, `X` if not a single *de novo* read is found and `Y:?` otherwise.
- `per_allele_reads_a` Number of reads partitioned per allele in sample A (allele1, allele2)
- `per_allele_reads_b` Number of reads partitioned per allele in sample B (allele1, allele2)
- `a_dropout` Coverage cut-off dropout detection using HP tags from phasing tools in sample A; possible values: Full dropout (`FD`), Haplotype dropout (`HD`), Not (`N`)
- `b_dropout` Coverage cut-off dropout detection using HP tags from phasing tools in sample B; possible values: Full dropout (`FD`), Haplotype dropout (`HD`), Not (`N`)
- `index` Index of this allele in the TRGT VCF
- `a_MC` TRGT VCF motif counts for this locus in sample A
- `b_MC` TRGT VCF motif counts for this locus in sample B
- `a_AL` TRGT VCF allele lengths for this locus in sample A
- `b_AL` TRGT VCF allele lengths for this locus in sample B
- `b_overlap_coverage` Reciprocal of `denovo_coverage`, the number of reads in per allele in sample B that overlap compared to the sample

## Individual read alignment scores

Both `trio` and `duo` accept `--alignment-scores <PATH>` to export the scores used during de novo assessment into a separate, headered TSV file:

```bash
trgt-denovo trio \
  --reference reference.fa --bed repeats.bed \
  --mother mother --father father --child child \
  --out calls.tsv --alignment-scores alignment_scores.tsv
```

This export does not change the main results table or the analysis. 

Each row represents one successfully completed read-to-candidate-allele alignment:

| Column | Meaning |
|---|---|
| `chrom` | Locus contig |
| `start` | Zero-based locus start, as in the BED and main results table |
| `end` | Exclusive locus end, as in the BED and main results table |
| `trid` | Tandem repeat identifier |
| `candidate_genotype` | Candidate allele's VCF allele index, matching `genotype` in the main results table |
| `candidate_index` | Candidate allele's zero-based genotype slot, matching `index` in the main results table |
| `sample_role` | `mother`, `father`, or `child` in trio mode; `a` or `b` in duo mode |
| `source_genotype` | VCF allele index of the allele to which the read was assigned |
| `source_index` | Zero-based genotype slot of the allele to which the read was assigned |
| `read_id` | Read name from the spanning BAM |
| `score` | Exact integer alignment score used by the analysis |

VCF allele indexes use `0` for the reference allele and `1`, `2`, etc. for alternate alleles. Genotype slots refer to the sample's genotype sorted by VCF allele index; they distinguish the two copies even for a homozygous genotype. For example, genotype `1/1` has two source slots, `0` and `1`, both with `source_genotype=1`. The candidate sample is always the child in trio mode or sample A in duo mode.

For each candidate child allele, trio mode exports scores for reads assigned to every maternal and paternal allele, plus the child reads assigned to that candidate allele. Duo mode exports scores for reads assigned to every B allele, plus the A reads assigned to that candidate allele. Scores are exported for all assessed candidates, including those ultimately classified as not de novo. A parental or B read can therefore appear in multiple rows, one for each candidate against which it was aligned.

Scores are non-positive penalties: `0` indicates no penalty in the scored alignment region, and higher (less negative) values indicate greater similarity.

Unsuccessful alignments have no row or placeholder score. Loci skipped before de novo alignment, including quick-mode skips and allele-loading failures such as missing genotypes, produce no score rows. If no alignments complete, the file still contains its header.
