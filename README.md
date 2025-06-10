<h1 align="center"></></h1>

<h1 align="center">TRGT-denovo</h1>

<h3 align="center">Calling <em>de novo</em> tandem repeat mutations</h3>

***

TRGT-denovo is a companion tool of [TRGT (Tandem Repeat Genotyper)](https://github.com/PacificBiosciences/trgt) that further annotates TRGT genotyping using additional read-level information from samples, specifically for identifying *de novo* tandem repeat mutations in parent-child trios and 1-to-1 comparisons using PacBio Hifi data.

Preprint: [TRGT-denovo: accurate detection of de novo tandem repeat mutations](https://doi.org/10.1101/2024.07.16.600745)

Developers: [Tom Mokveld](https://github.com/tmokveld), [Egor Dolzhenko](https://github.com/egor-dolzhenko)

## Early version warning

Please note that TRGT-denovo is still in development and is subject to significant changes that can affect anything from input / output file formats to program behavior.

## Availability

* [Latest release with binary](https://github.com/PacificBiosciences/trgt-denovo/releases/latest)

To build TRGT-denovo you need a working C compiler: it was tested on Linux with Clang 13.0.0 & GCC 11.3.0 and on Mac OSX (M1) with Clang 15.0.7 & GCC 14.0.0.

## Documentation

* [Getting started](docs/example.md)
* [Output format](docs/output.md)
* [Interpretation](docs/interpretation.md)
* [Command-line interface](docs/cli.md)

## Need help?
If you notice any missing features, bugs, or need assistance with analyzing the output of TRGT-denovo, 
please don't hesitate to open a GitHub issue.

## Support information
TRGT-denovo is software intended for research use only and not for use in diagnostic procedures. 
While efforts have been made to ensure that TRGT-denovo lives up to the quality that PacBio strives for, we make no warranty regarding this software.

As TRGT-denovo is not covered by any service level agreement or the like, please do not contact a PacBio Field Applications Scientists or PacBio Customer Service for assistance with any TRGT-denovo release. 
Please report all issues through GitHub instead. 
We make no warranty that any such issue will be addressed, to any extent or within any time frame.

## Changelog

- 0.2.3
  - Bug fix: No longer skip sites that are actually genotyped when tandem repeats overlap

- 0.2.2
  - Output now includes the columns chrom, start, end, and motifs.
  - Identical starting positions are now allowed for tandem repeats, with reads pruned if they do not match the locus TRID.
  - Read IDs contributing to de novo coverage are now logged, enabling traceability of TRIDs and alleles. The log file is named based on trio or duo mode: {child_name}_trio_denovo_reads.txt or {a_name}_duo_denovo_reads.txt.
  - Bug fix: Resolved a flag conflict in duo mode by adjusting short-form flags, preventing overlap between --bed (-b) and --sample-b.

- 0.2.1
  - Duo and trio mode now always output entries regardless of error status (e.g., missing genotyping, skip because of quick mode).
  - Add `denovo_status` analog to duo mode.

- 0.2.0
  - Implemented duo mode, it is now possible to perform 1-to-1 sample comparisons, following the same principles as in trio analysis. This can be done using the subcommand `trgt-denovo duo`.
  - Implemented the `--quick` flag. Users can now specify `--quick AL[,<fraction>]` to skip loci where allele lengths are similar between parents and child or between two samples. If no fraction is specified (or fraction is 0), it checks for exact matches. If a fraction is specified, it checks if the relative difference is within the given tolerance.
  - Added a Jupyter notebook in scripts/python/trio_analysis.ipynb to describe a simple trio analysis to do *de novo* candidate selection.

- 0.1.3
  - Changes to TRGT-denovo output:
      - Truncate zeros in output.
      - Report TRGT allele lengths observed in each family member as `sample_AL`.
      - Renamed TRGT motif counts from `sample_motif_counts` to `sample_MC`.
      - Report the overlap coverage; this is the reciprocal of the de novo coverage, i.e., the number of reads in per allele in the parent that overlap compared to the child data.
  - Lower memory footprint: Better memory management, significantly reduces memory usage with large repeat catalogs.
  - Improved IO error handling.

- 0.1.2
  - With the recent changes to TRGT, within-sample partitioning is now alignment-free in TRGT-denovo. An optional parameter has been added to still use alignment (`--partition-by-aln`).
  - Homozygous alleles are no longer collapsed: de novo evidence will now always gathered and be specific to a single allele only.

- 0.1.1
  - Add cli parameter to set aligner penalties.
  - Document the codebase.
  - Update documentation to include interpretation of generated output.

## Ongoing work

- Add downstream scripts to: 
  - further annotate TRGT-denovo output with population specific data, e.g., outlier detection, is a particular allele an outlier with respect to a population.
  - select candidate de novo calls through classification.
  - show an example of duo based analysis to detect sites of interest.
- Haplotype matching across samples using flanking variation.

### DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.