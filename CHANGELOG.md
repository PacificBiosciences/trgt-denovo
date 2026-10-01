# Changelog

Any changes to TRGT-denovo are noted below:

## [0.4.0]

- Add `--alignment-scores <PATH>` to `trio` and `duo` to export per-read alignment scores as a separate TSV.
- Add `scripts/python/trio_alignment_scores.ipynb` to visualize trio alignment scores.

## [0.3.0]

- Switched file handling to `rust-htslib`.
- CLI now also accepts explicit BAM/VCF file paths for each sample (`--mother-vcf/--mother-bam`, `--sample-a-vcf/--sample-a-bam`, etc.), alongside prefix-based inputs.
- Bug fix: normalize observed genotypes such that the phasing state does not affect allele sequence ordering, leading to unexpected results.

## [0.2.3]

- Bug fix: No longer skip sites that are actually genotyped when tandem repeats overlap.

## [0.2.2]

- Output now includes the columns chrom, start, end, and motifs.
- Identical starting positions are now allowed for tandem repeats, with reads pruned if they do not match the locus TRID.
- Read IDs contributing to de novo coverage are now logged, enabling traceability of TRIDs and alleles. The log file is named based on trio or duo mode: {child_name}_trio_denovo_reads.txt or {a_name}_duo_denovo_reads.txt.
- Bug fix: Resolved a flag conflict in duo mode by adjusting short-form flags, preventing overlap between --bed (-b) and --sample-b.

## [0.2.1]

- Duo and trio mode now always output entries regardless of error status (e.g., missing genotyping, skip because of quick mode).
- Add `denovo_status` analog to duo mode.

## [0.2.0]

- Implemented duo mode, it is now possible to perform 1-to-1 sample comparisons, following the same principles as in trio analysis. This can be done using the subcommand `trgt-denovo duo`.
- Implemented the `--quick` flag. Users can now specify `--quick AL[,<fraction>]` to skip loci where allele lengths are similar between parents and child or between two samples. If no fraction is specified (or fraction is 0), it checks for exact matches. If a fraction is specified, it checks if the relative difference is within the given tolerance.
- Added a Jupyter notebook in scripts/python/trio_analysis.ipynb to describe a simple trio analysis to do *de novo* candidate selection.

## [0.1.3]

- Changes to TRGT-denovo output:
    - Truncate zeros in output.
    - Report TRGT allele lengths observed in each family member as `sample_AL`.
    - Renamed TRGT motif counts from `sample_motif_counts` to `sample_MC`.
    - Report the overlap coverage; this is the reciprocal of the de novo coverage, i.e., the number of reads in per allele in the parent that overlap compared to the child data.
- Lower memory footprint: Better memory management, significantly reduces memory usage with large repeat catalogs.
- Improved IO error handling.

## [0.1.2]

- With the recent changes to TRGT, within-sample partitioning is now alignment-free in TRGT-denovo. An optional parameter has been added to still use alignment (`--partition-by-aln`).
- Homozygous alleles are no longer collapsed: de novo evidence will now always gathered and be specific to a single allele only.

## [0.1.1]

- Add cli parameter to set aligner penalties.
- Document the codebase.
- Update documentation to include interpretation of generated output.
