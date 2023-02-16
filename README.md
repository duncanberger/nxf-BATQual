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

- __fastq_files.csv__: Comma delineated list of fastq files, in the format: 'sample_id,read1,read2' (with this header included)
  - `sample`: sample ID
  - `read1`: forward reads
  - `read2`: reverse reads
- __fasta_files.csv__: Comma delineated list of fastq files, in the format: sample_id,fasta (with this header included)
  - `sample`: sample ID
  - `fasta`: fasta files

- __Genome sequence__: FASTA file containing the genome sequence. Can be retrieved from NCBI.
- __Gene annotation file__: GFF file containing the genome annotation. Can be retrieved from NCBI.
- __Sample file__: TSV file containing sample information. Must contain the following columns:
  - `sample`: sample ID
  - `file_name`: name of the FASTQ file.
  - `group`: grouping factor for differential expression and exploratory plots.
  - `rep_no`: repeat number (if more than one sample per group).

Explanation of parameters:
- `ref_genome`: genome sequence for mapping reads.
- `ref_ann`: annotation of genes/features in the reference genome.
- `sample_file`: TSV file containing sample information (see below)
- `data_dir`: path to directory containing FASTQ files.
- `paired`: data are paired-end (default is to assume single-end)
- `strandedness`: is data stranded? Options: `unstranded`, `forward`, `reverse`. Default = `reverse`.
- `cont_tabl`: (optional) table of contrasts to be performed for differential expression.
- `func_file`: (optional) functional annotation file - if provided, functional enrichment of DE genes will be performed.
- `p_thresh`: adjusted p-value threshold for identifying differentially expressed genes. Default = 0.05.
- `l2fc_thresh`: absolute log2(FoldChange) threshold for identifying differentially expressed genes. Default = 1.
- `skip_trimming`: do not trim adaptors from reads.
- `outdir`: the output directory where the results will be saved (Default: `./results`).




  Example:

  If data are single-end, leave the `file2` column blank.

    ```console
    sample	file1   file2	group	rep_no
    AS_1	SRX1607051_T1.fastq.gz	    Artificial_Sputum	1
    AS_2	SRX1607052_T1.fastq.gz	    Artificial_Sputum	2
    AS_3	SRX1607053_T1.fastq.gz	    Artificial_Sputum	3
    MB_1	SRX1607054_T1.fastq.gz	    Middlebrook	1
    MB_2	SRX1607055_T1.fastq.gz	    Middlebrook	2
    MB_3	SRX1607056_T1.fastq.gz	    Middlebrook	3
    ER_1	SRX1607060_T1.fastq.gz	    Erythromycin	1
    ER_2	SRX1607061_T1.fastq.gz	    Erythromycin	2
    ER_3	SRX1607062_T1.fastq.gz	    Erythromycin	3
    KN_1	SRX1607066_T1.fastq.gz	    Kanamycin	1
    KN_2	SRX1607067_T1.fastq.gz	    Kanamycin	2
    KN_3	SRX1607068_T1.fastq.gz	    Kanamycin	3
    ```

## Output

1. __trim_galore__ directory containing adaptor-trimmed RNA-Seq files and FastQC results.
2. __read_counts__ directory containing:
    1. `ref_gene_df.tsv`: table of genes in the annotation.
    2. `gene_counts.tsv`: raw read counts per gene.
    3. `cpm_counts.tsv`: size factor scaled counts per million (CPM).
    4. `rpkm_counts.tsv`: size factor scaled and gene length-scaled counts, expressed as reads per kilobase per million mapped reads (RPKM).
3. __PCA_samples__ directory containing principal component analysis results.
4. __diff_expr__ directory containing differential expression results.
5. __func_enrich__ directory containing functional enrichment results (optional).

