# Interpreting TRGT-denovo output

Trio fields:

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
- `denovo_status` Indicates if the allele is *de novo*, only if `allele_origin` is defined; possible values: `{X, Y:{+, -, =}}`. This is `X` if no *de novo* read is found and `Y` otherwise, if parental origin can be determined without ambiguity the allele sequences can be compared directly such that the *de novo* type can be established as `+` (expansion), `-` (contraction), or `=` (substitution)
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

- `trid` ID of the tandem repeat, encoded as in the BED file
- `genotype` Genotype ID of sample A for a specific allele, corresponds to the TRGT genotype ID
- `denovo_coverage` Number of sample A reads supporting a *de novo* allele compared to sample B
- `allele_coverage` Number of sample A reads mapped to this specific allele
- `allele_ratio` Ratio of *de novo* coverage to allele coverage
- `a_coverage` Total number of sample A reads at this site
- `a_ratio` Ratio of *de novo* coverage to total coverage at this site
- `mean_diff_b` Score difference between *de novo* and sample B reads; lower values indicate greater similarity
- `per_allele_reads_a` Number of reads partitioned per allele in sample A (allele1, allele2)
- `per_allele_reads_b` Number of reads partitioned per allele in sample B (allele1, allele2)
- `a_dropout` Coverage cut-off dropout detection using HP tags from phasing tools in sample A; possible values: Full dropout (`FD`), Haplotype dropout (`HD`), Not (`N`) 
- `b_dropout` Coverage cut-off dropout detection using HP tags from phasing tools in sample B; possible values: Full dropout (`FD`), Haplotype dropout (`HD`), Not (`N`) 
- `index` Index of this allele in the TRGT VCF
- `a_MC` TRGT VCF motif counts for this locus in sample A
- `b_MC` TRGT VCF motif counts for this locus in sample B
- `a_AL` TRGT VCF allele lengths for this locus in sample A
- `b_AL` TRGT VCF allele lengths for this locus in sample B
- `b_overlap_coverage` Reciprocal of `denovo_coverage`, the number of reads in per allele in sample B that overlap compared to the sample A data