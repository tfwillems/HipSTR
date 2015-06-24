# HipSTR
**H**aplotype-based **i**mputation, **p**hasing and genotyping of **STR**s

#### Author: Thomas Willems <twillems@mit.edu>
#### License: GNU v2

## Introduction
Short tandem repeats [(STRs)](http://en.wikipedia.org/wiki/Microsatellite) are highly repetitive genomic sequences comprised of repeated copies of an underlying motif. Prevalent in most organisms' genomes, STRs are of particular interest because they mutate much more rapidly than most other genomic elements. As a result, they're extremely informative for genomic identification, ancestry inference and genealogy.

Despite their utility, STRs are particularly difficult to genotype . The repetitive sequence responsible for their high mutability also results in frequent alignment errors that can complicate and bias downstream analyses. In addition, PCR stutter errors often result in reads that contain additional or fewer repeat copies than the true underlying genotype. 

**HipSTR** was specifically developed to deal with these errors in the hopes of obtaining more robust STR genotypes. In particular, it accomplishes this by:

1. Learning locus-specific PCR stutter models using an [EM algorithm] (http://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm)
2. Mining candidate STR alleles from population-scale sequencing data
2. Utilizing phased SNP haplotypes to genotype, phase and/or impute STRs
3. Employing a specialized hidden Markov model to align reads to candidate sequences while accounting for stutter



## Installation
HipSTR requires a standard c++ compiler as well as Java version 1.7 or later.
To obtain HipSTR and all of its associated  submodules, use:

    % git clone --recursive https://github.com/tfwillems/HipSTR.git

To build, use Make:

    % cd vcflib
    % make

On Mac, before running Make, change the line in *vcflib/smithwaterman/Makefile* from

    % LDFLAGS=-Wl,-s
to

    % LDFLAGS=-Wl


## Usage
### De novo STR calling
    % ./HipSTR --bams    sample1.bam,sample2.bam,sample3.bam 
               --indexes sample1.bam.bai,sample2.bam.bai,sample3.bam.bai
               --fasta   /data/hg19_by_chrom/
               --regions /data/str_regions.bed
               --str-vcf str_calls.vcf.gz

### STR calling using a reference panel of STRs


### STR imputation




### Alignment Visualization
