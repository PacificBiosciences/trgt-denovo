# Interpreting TRGT-denovo output

TRGT-denovo scores and then outputs target sites in a tab-separated format, with 11 columns:

- `trid` ID of the tandem repeat, encoded as in the BED file
- `genotype` Genotype ID of the child for a specific allele, corresponds to the TRGT genotype ID
- `denovo_coverage` Number of child reads supporting a de novo allele compared to parental data.
- `allele_coverage` Number of child reads mapped to this specific allele
- `allele_ratio` Ratio of de novo coverage to allele coverage.
- `child_coverage` Total number of child reads at this site
- `child_ratio` Ratio of de novo coverage to total coverage at this site
- `mean_diff_father` Score difference between de novo and paternal reads; lower values indicate greater similarity
- `mean_diff_mother` Score difference between de novo and maternal reads; lower values indicate greater similarity
- `father_dropout_prob` Dropout rate for reads coming from the mother
- `mother_dropout_prob` Dropout rate for reads coming from the father.
- `allele_origin` Inferred origin of the allele based on alignment; possible values: `{F:{1,2,?}, M:{1,2,?}, ?}`
- `denovo_status` Indicates if the allele is de novo, only if `allele_origin` is defined; possible values: `{?, Y:{+, -, =}}`
- `per_allele_reads_father` Number of reads partitioned per allele in the father (allele1, allele2)
- `per_allele_reads_mother` Number of reads partitioned per allele in the mother (allele1, allele2)
- `per_allele_reads_child` Number of reads partitioned per allele in the child (allele1, allele2)
- `father_dropout` Coverage cut-off dropout detection in father; possible values: Full dropout (`FD`), Haplotype dropout (`HD`), Not (`N`) 
- `mother_dropout` Coverage cut-off dropout detection in mother; possible values: Full dropout (`FD`), Haplotype dropout (`HD`), Not (`N`) 
- `child_dropout` Coverage cut-off dropout detection in child; possible values: Full dropout (`FD`), Haplotype dropout (`HD`), Not (`N`) 
- `index` Index of this allele in the TRGT VCF, used for linking to `child_motif_counts`
- `father_motif_counts` TRGT VCF motif counts for this locus in the father
- `mother_motif_counts` TRGT VCF motif counts for this locus in the mother
- `child_motif_counts` TRGT VCF motif counts for this locus in the child
