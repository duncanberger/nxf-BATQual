# nxf-bact_typ


[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)


## Introduction

**nxf-bact_typ** is a Nextflow pipeline for performing assembling, annotation and typing bacterial genomes. 

## Pipeline summary

When processing Illumina sequencing reads ('--mode fastq'), the pipeline will perform the following steps:

1. Trim adaptors from reads ([`fastp`](https://github.com/OpenGene/fastp)).
2. Assemble one genome per-isolate using ([`VelvetOptimiser`](https://github.com/tseemann/VelvetOptimiser)).
3. Estimate heterozygosity with the read data to check for contamination using ([`minimap`](https://github.com/lh3/minimap2)) and ([`BCFtools`](https://samtools.github.io/bcftools/bcftools.html)).
4. Calculate standard assembly statistics
5. Estimate assembly completeness using ([`BUSCO`](https://busco.ezlab.org/))
6. Check for assembly completeness and contamination using ([`CheckM`](https://github.com/Ecogenomics/CheckM))
7. Estimate most likely species based on the closet ([`MASH`](https://github.com/marbl/Mash)) hits in the RefSeq database. 
8. Perform read-based pneumococcal serotyping using ([`seroBA`](https://github.com/sanger-pathogens/seroba))
9. Perform assembly-based pneumococcal serotyping using ([`pneumoKITy`](https://github.com/sanger-pathogens/seroba))
10. Assign Global Pneumococcal Sequence Clusters ([`GPSCs`](https://www.pneumogen.net/gps/)) using popunk ([`GPSCs`](https://poppunk.net/)) 
11. Annotate each genome using ([`Prokka`](https://github.com/tseemann/prokka)) 

When processing genome assemblies ('--mode fasta'), the pipeline will omit read-specific steps:
1. Calculate standard assembly statistics
2. Estimate assembly completeness using ([`BUSCO`](https://busco.ezlab.org/))
3. Check for assembly completeness and contamination using ([`CheckM`](https://github.com/Ecogenomics/CheckM))
4. Estimate most likely species based on the closet ([`MASH`](https://github.com/marbl/Mash)) hits in the RefSeq database. 
5. Perform assembly-based pneumococcal serotyping using ([`pneumoKITy`](https://github.com/sanger-pathogens/seroba))
6. Assign Global Pneumococcal Sequence Clusters ([`GPSCs`](https://www.pneumogen.net/gps/)) using popunk ([`GPSCs`](https://poppunk.net/)) 
7. Annotate each genome using ([`Prokka`](https://github.com/tseemann/prokka)) 

## Installation

You will need to install [`Nextflow`](https://www.nextflow.io/) (version 21.10.3+).

You can run the pipeline as follows:

    nextflow run main.nf --input fastq_files.csv --mode fastq
or 
    nextflow run main.nf --input fasta_files.csv --mode fasta

The `-resume` parameter will re-start the pipeline if it has been previously run.

## Required input

- __fastq_files.csv__: Comma delineated list of fastq files, in the format: 'sample_id,read1,read2'. A header 'sample_id,read1,read2' must be included. 
  - `sample`: sample ID
  - `read1`: forward reads
  - `read2`: reverse reads
- __fasta_files.csv__: Comma delineated list of fastq files, in the format: 'sample_id,fasta'. A header 'sample_id,fasta' must be included. 
  - `sample`: sample ID
  - `fasta`: fasta files

## Optional input

### Velvet assembly
`--min_k` : Minimum k-mer length for velvet genome assembly [90]. <br />
`--max_k` : Minimum k-mer length for velvet genome assembly [141].

### BUSCO 
`--lineage` : Specify the BUSCO lineage to be used ["lactobacillales_odb10"].

### Assembly filtering
`--minContigLength` : Filter for minimum contig length in output [0].

### Output
`--outdir` : The output directory where the results will be saved ["results"]. <br />
`--name` : Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

## Pipeline results

#1. __trim_galore__ directory containing adaptor-trimmed RNA-Seq files and FastQC results.
#2. __read_counts__ directory containing:
#    1. `ref_gene_df.tsv`: table of genes in the annotation.
#    2. `gene_counts.tsv`: raw read counts per gene.
#    3. `cpm_counts.tsv`: size factor scaled counts per million (CPM).
#    4. `rpkm_counts.tsv`: size factor scaled and gene length-scaled counts, expressed as reads per kilobase per million mapped reads (RPKM).
#3. __PCA_samples__ directory containing principal component analysis results.
#4. __diff_expr__ directory containing differential expression results.
#5. __func_enrich__ directory containing functional enrichment results (optional).

