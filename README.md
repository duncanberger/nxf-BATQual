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


