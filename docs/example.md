# Introduction

A brief overview of the steps needed to call *de novo* tandem repeat mutations using TRGT-denovo given TRGT output.

## Prerequisites

- [TRGT binary](https://github.com/PacificBiosciences/trgt/releases/latest)
- [TRGT-denovo binary](https://github.com/PacificBiosciences/trgt-denovo/releases/latest)
- Repeat definition files are available in [this Zenodo repository](https://zenodo.org/record/8329210)
  and definitions of known pathogenic repeats are [also available here](https://github.com/PacificBiosciences/trgt/tree/main/repeats/) (always the same sites or a subset of those used with TRGT)
- Aligned HiFi data of a family trio (father, mother, and child)
- The reference genome used for read alignment




## Calling *de novo* tandem repeat mutations

Given the following data: 

- Reference genome `reference.fasta`
- Repeat definition file `repeat.bed`.
- Aligned sequencing data of the family (father, mother, and son respectively) `sample_F.bam`, `sample_M.bam`, and `sample_S.bam`.

### Data pre-processing

All data must first be genotyped by TRGT:

```
./trgt --genome reference.fasta \
       --repeats repeat.bed \
       --reads sample_F.bam \
       --output-prefix sample_F \
       --karyotype XY
```

```
./trgt --genome reference.fasta \
       --repeats repeat.bed \
       --reads sample_M.bam \
       --output-prefix sample_M \
       --karyotype XX
```

```
./trgt --genome reference.fasta \
       --repeats repeat.bed \
       --reads sample_S.bam \
       --output-prefix sample_S \ 
       --karyotype XY
```

TRGT outputs the genotyped repeat sites in a VCF file stored in `prefix.vcf.gz` and the spanning reads that were used to genotype each site (that fully span the repeat sequences) stored in `prefix.spanning.bam`. TRGT-denovo requires sorted BAM and VCF data, hence you will need to sort and index the output VCF and BAM files. For each family member this involves:

#### VCF sorting
```
bcftools sort -Ob -o sample_F.sorted.vcf.gz sample_F.vcf.gz
bcftools index sample_F.sorted.vcf.gz
```

#### BAM sorting
```
samtools sort -o sample_F.spanning.sorted.bam sample_F.spanning.bam
samtools index sample_F.spanning.sorted.bam
```

Such that you end up with `sample_F.sorted.vcf.gz`, `sample_F.spanning.sorted.bam`, `sample_M.sorted.vcf.gz`, `sample_M.spanning.sorted.bam`, `sample_S.sorted.vcf.gz`, `sample_S.spanning.sorted.bam` (and their associated `.bam.bai` and `.vcf.gz.csi` indices).

### Running TRGT-denovo

With all preprocessing completed, we can call *de novo* repeat expansion mutations using TRGT-denovo from the sample data. Note that family members are supplied by their common prefix of `spanning.sorted.bam` and `sorted.vcf.gz`, i.e., `sample_F`, `sample_M`, and `sample_S` and path if not running TRGT-denovo in the same directory as the data:

```
./TRGT-denovo trio --reference reference.fasta \
              --bed repeat.bed \
              --father sample_F \
              --mother sample_M \
              --child sample_S \
              --out out.csv
```

For further interpretation of TRGT-denovo output see [here](interpretation.md)