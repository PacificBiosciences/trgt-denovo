# Introduction

A brief overview of the steps necessary to perform *de novo* tandem repeat mutation calling from TRGT output using TRGT-denovo

## Prerequisites

- [TRGT binary](https://github.com/PacificBiosciences/trgt/releases/latest)
- [TRGT-denovo binary](https://github.com/PacificBiosciences/trgt-denovo/releases/latest)
- Aligned HiFi data of a family trio (father, mother, child)
- The reference genome used for the alignments
- A BED file with repeat expansion definitions (always same or a subset of those with TRGT)

## Calling *de novo* tandem repeats

Given the following data: 

- reference genome `reference.fasta`
- aligned sequencing data of the family (father, mother, and son respectively) `sample_F.bam`, `sample_M.bam`, `sample_S.bam`, 
- repeat definition file `repeat.bed`.

### Data pre-processing

Prior to *de novo* calling all data must first be genotyped by TRGT:

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

With all preprocessing completed it we can call *de novo* repeat expansion mutations using TRGT-denovo from the sample data. Note that family members are supplied by their common prefix of `spanning.sorted.bam` and `sorted.vcf.gz`, i.e., `sample_F`, `sample_M`, and `sample_S` and path if not running TRGT-denovo in the same directory as the data:

```
./TRGT-denovo trio --reference reference.fasta \
              --bed repeat.bed \
              --father sample_F \
              --mother sample_M \
              --child sample_S \
              --out out.csv
```

## Example

Below output of HG002 is shown for two candidate *de novo* tandem repeat mutation sites:
```
trid	genotype	denovo_coverage	allele_coverage	allele_ratio	child_coverage	child_ratio	mean_diff_father	mean_diff_mother	father_dropout_prob	mother_dropout_prob	allele_origin	denovo_status	per_allele_reads_father	per_allele_reads_mother	per_allele_reads_child	index	father_motif_counts	mother_motif_counts	child_motif_counts	maxlh
chr1_47268728_47268830_ATAA	1	19	37	0.5135	37	0.5135	6.7368	6.7368	0.0000	0.0000	M:1	Y:=	43	26	37	0	25	25	25	0.7152
chr1_7862944_7863157_TATTG	1	0	21	0.0000	37	0.0000	0.0000	19.2000	0.0000	0.0000	F:2	X	18,17	16,19	21,16	0	27,29	27,63	29,60	0.9820
chr1_7862944_7863157_TATTG	2	16	16	1.0000	37	0.4324	171.8750	22.8750	0.0000	0.0000	M:2	Y:-	18,17	16,19	21,16	1	27,29	27,63	29,60	1.0000
```

### Site 1

The first site is homozygous in the child, hence only one call is made. It has a *de novo* coverage (DNC) of 19, i.e., there are 19 reads that support a candidate *de novo* allele relative to the parental read alignments. The DNC should always be put into context of the total coverage, to ascertain that:

1. There is sufficient coverage
2. The ratio between the two is close to 0.5

The DNC at this site is high and the ratio makes it likely that this is a confident call. However, the score difference with respect to either parents is low, such that the expected event size is small. Additionally, the score difference is equivalent in both parents such that parental origin may not be assessed. Generally the parent with the smallest score difference is the one from which is inherited:

![Example site 1](figures/example_site_1.png)

### Site 2

The second site is heterozygous in the child, hence both child alleles are tested. The second allele is a potential *de novo* call. The score difference with respect to the maternal alleles is the smallest.

![Example site 2](figures/example_site_2.png)

