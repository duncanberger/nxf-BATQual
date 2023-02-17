# nxf-BATQual [Pre-release]
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.4-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)

## Table of contents 
* [Introduction](#Introduction)
* [Pipeline summary](#pipeline_summary)
* [Installation](#install)
* [Running nxf-BATQual](#run)
* [Aggregation and quality control](#aggregate)


## Introduction <a name="Introduction"></a>

**[NOTE: CURRENTLY IN ACTIVE DEVELOPMENT]** <br />
**nxf-BATQual** is a Nextflow pipeline for performing assembling, annotation and typing bacterial genomes.

## Pipeline summary <a name="pipeline_summary"></a>
When processing Illumina sequencing reads ('--mode fastq'), the pipeline will perform the following steps:

1. Trim adaptors from reads ([`fastp`](https://github.com/OpenGene/fastp)).
2. Assemble one genome per-isolate using ([`VelvetOptimiser`](https://github.com/tseemann/VelvetOptimiser)).
3. Estimate heterozygosity with the read data to check for contamination using ([`minimap`](https://github.com/lh3/minimap2)) and ([`BCFtools`](https://samtools.github.io/bcftools/bcftools.html)).
4. Calculate standard assembly statistics
5. Estimate assembly completeness using ([`BUSCO`](https://busco.ezlab.org/))
6. Check for assembly completeness and contamination using ([`CheckM`](https://github.com/Ecogenomics/CheckM))
7. Evaluate assembly using ([`QUAST`](https://quast.sourceforge.net/)) 
8. Estimate most likely species based on the closet ([`MASH`](https://github.com/marbl/Mash)) hits in the RefSeq database. 
9. Perform read-based pneumococcal serotyping using ([`seroBA`](https://github.com/sanger-pathogens/seroba))
10. Perform assembly-based pneumococcal serotyping using ([`pneumoKITy`](https://github.com/sanger-pathogens/seroba))
11. Assign Global Pneumococcal Sequence Clusters ([`GPSCs`](https://www.pneumogen.net/gps/)) using popunk ([`PopPunk`](https://poppunk.net/)) 
12. Perform multilocus sequence typing using ([`mlst`](https://github.com/tseemann/mlst))
13. Annotate each genome using ([`Prokka`](https://github.com/tseemann/prokka)) 

When processing genome assemblies ('--mode fasta'), the pipeline will omit read-specific steps:
1. Calculate standard assembly statistics
2. Estimate assembly completeness using ([`BUSCO`](https://busco.ezlab.org/))
3. Check for assembly completeness and contamination using ([`CheckM`](https://github.com/Ecogenomics/CheckM))
4. Evaluate assembly using ([`QUAST`](https://quast.sourceforge.net/)) 
5. Estimate most likely species based on the closet ([`MASH`](https://github.com/marbl/Mash)) hits in the RefSeq database. 
6. Perform assembly-based pneumococcal serotyping using ([`pneumoKITy`](https://github.com/sanger-pathogens/seroba))
7. Assign Global Pneumococcal Sequence Clusters ([`GPSCs`](https://www.pneumogen.net/gps/)) using popunk ([`PopPunk`](https://poppunk.net/)) 
8. Perform multilocus sequence typing using ([`mlst_check`](https://github.com/sanger-pathogens/mlst_check))
9. Annotate each genome using ([`Prokka`](https://github.com/tseemann/prokka)) 

### Schematic overview 
![rect13856](https://user-images.githubusercontent.com/29282405/219794605-086076a0-2e6d-471b-8e48-501d00f853b2.png)

## Installation <a name="install"></a>
### Software
You will need to install [`Nextflow`](https://www.nextflow.io/) (version 21.10.3+).
You will need to install [`pneumoKITy`](https://github.com/CarmenSheppard/PneumoKITy) manually.
```
# Clone this repo
git clone https://github.com/duncanberger/nxf-BATQual.git

# Set up conda dependencies 
conda create -n nxf-BATQual -f environment.yml 
```
### Databases
This pipeline uses a number of databases, before running the pipeline you should check that version included below are the most relevant/up to date for your analyses.  The can be downloaded and made ready for analysis as follow:

#### RefSeq Mash database
```
# Download a preformatted MASH database of RefSeq (release 70, if you want a more up to date version you'll need to build your own)
wget https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh

# Move it to the database folder
mv refseq.genomes.k21s1000.msh DB/
```
#### seroBA
```
# seroBA can be found here: https://github.com/sanger-pathogens/seroba.git but all we want is the database folder. You can either clone the whole repo
git clone https://github.com/sanger-pathogens/seroba.git

# Or you can download that folder specifically
git init
git remote add origin -f https://github.com/sanger-pathogens/seroba.git
echo "database" > .git/info/sparse-checkout
git config core.sparseCheckout true
git pull origin master

# In either case just move the database folder into the DB/ folder
mv database DB/
```
#### Global Pneumococcal Sequencing (GPS) clusters
```
# GPS reference database (n=42,163):
wget https://gps-project.cog.sanger.ac.uk/GPS_v6.zip

# GPS designation (933 GPSCs):
wget https://www.pneumogen.net/gps/GPS_v6_external_clusters.csv

mv GPS_v6* DB/
cd DB/
unzip GPS_v6.zip
```
#### Reference genome for the target species
```
# Download any relevant reference assembly for your species. E.g: 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/026/665/GCA_000026665.1_ASM2666v1/GCA_000026665.1_ASM2666v1_genomic.fna.gz 
gunzip GCA_000026665.1_ASM2666v1_genomic.fna.gz
mv GCA_000026665.1_ASM2666v1_genomic.fna DB/

# Then add the assembly name to 'main.config' file at the 'ref_assembly=' row
```

#### Benchmarking Universal Single-Copy Orthologs (BUSCO)
```
# By default BUSCO will download the relevant databases when trying to run, however, when running in nextflow it'll try to do that for every sample. A better solution is to download the relevant database prior to running the pipeline and specifying the relevant path. 
wget https://busco-data.ezlab.org/v5/data/lineages/lactobacillales_odb10.2020-03-06.tar.gz
mv lactobacillales_odb10.2020-03-06.tar.gz DB/
tar xvf lactobacillales_odb10.2020-03-06.tar.gz
```
#### Kraken 2 database
```
# Download the Kraken DB (capped at 8 Gb):
https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20221209.tar.gz

# Move to database directory and unpack
mv k2_standard_08gb_20221209.tar.gz DB/
cd DB/
tar xvf k2_standard_08gb_20221209.tar.gz
```

## Running nxf-BATQual <a name="run"></a>

You can run the pipeline as follows:
```
nextflow run main.nf --input fastq_files.csv --mode fastq
```
or 
```
nextflow run main.nf --input fasta_files.csv --mode fasta
```

The `-resume` parameter will re-start the pipeline if it has been previously run.


### Required input
- __fastq_files.csv__: Comma delineated list of fastq files, in the format: 'sample_id,read1,read2'. A header 'sample_id,read1,read2' must be included. 
  - `sample`: sample ID
  - `read1`: forward reads
  - `read2`: reverse reads
- __fasta_files.csv__: Comma delineated list of fastq files, in the format: 'sample_id,fasta'. A header 'sample_id,fasta' must be included. 
  - `sample`: sample ID
  - `fasta`: fasta files

### Optional input
#### Velvet assembly
- `--min_k` : Minimum k-mer length for velvet genome assembly [90]. <br />
- `--max_k` : Minimum k-mer length for velvet genome assembly [141].

#### BUSCO 
- `--lineage` : Specify the BUSCO lineage to be used ["lactobacillales_odb10"].

#### Assembly filtering
- `--minContigLength` : Filter for minimum contig length in output [0].

#### Output
- `--outdir` : The output directory where the results will be saved ["results"]. <br />
- `--name` : Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

### Pipeline results
The output will be written to results/$sample_ID/* which will contain the following results files:

1. `$sample_ID.velvet_contigs.fa` : VelvetOptimiser genome assembly.
2. `$sample_ID.prokka.fna` : VelvetOptimiser genome assembly renamed by Prokka to match GFF file. 
3. `$sample_ID.prokka.gff` : VelvetOptimiser genome annotation, in Prokka format.
4. `$sample_ID.prokka_results.txt` : Counts of predicted genes. 
5. `$sample_ID.hetperc_results.txt` : Estimated heterozygosity of sequencing reads.
6. `$sample_ID.assembly_results.txt` : Assembly statistics.
7. `$sample_ID.BUSCO_results.txt` : seroBA serotyping results 
8. `$sample_ID.checkM_results.txt` : seroBA serotyping results .
9. `$sample_ID.mash_results.txt` : seroBA serotyping results .
10. `$sample_ID.seroBA_results.txt` : seroBA serotyping results .
11. `$sample_ID.pneumoKITy_results.txt` : seroBA serotyping results .
12. `$sample_ID.GPSC_results.txt` : seroBA serotyping results 

Each file is formatted into columns as follows: 
1. `sample_ID`: Isolate ID
2. `metric`: Metric used for evaluation
3. `results`: Metric results
4. `status`: QC status, unpopulated except for pneumoKITy QC, to be populated at the next stage. 

## Aggregation and quality control <a name="aggregate"></a>
To aggregate the results across multiple isolates and runs, I have written an accessory script to produce merged output tables and plot the results.

You can run the script as follows (where results is the name of the folder specified by the '--output' parameter:
```
python scripts/filter.py --input results
```
### Optional input
#### Output parameter
- `--output` : Output file name [aggregated_stats.* ]

#### QC parameters
I've added extensive QC parameters, genomes not meeting these parameters will be marked as 'FAIL' in the aggregated summaries.

- `--completeness_threshold` :            CheckM completeness threshold [99] <br />
- `--contamination_threshold` :           CheckM contamination threshold [1] <br />
- `--strain_heterogeneity_threshold` :    CheckM heterogeneity threshold [0.1] <br />
- `--busco_completeness_threshold` :      BUSCO completeness threshold (%) [95] <br />
- `--busco_duplication_threshold` :       BUSCO duplication threshold (%) [1] <br />
- `--busco_fragmented_threshold` :        BUSCO fragmented threshold (%) [1] <br />
- `--busco_missing_threshold` :           BUSCO missing threshold (%) [1] <br />
- `--assembly_length_threshold_min` :     Assembly length, lower threshold (bp) [1945246] <br />
- `--assembly_length_threshold_max` :     Assembly length, upper threshold (bp) [2255392] <br />
- `--gc_threshold_min` :                  GC content, lower threshold (%) [39.2] <br />
- `--gc_threshold_max` :                  GC content, upper threshold (%) [40] <br />
- `--gap_sum_threshold` :                 Total gap length (bp) [5631] <br />
- `--gap_count_threshold` :               Total gap count [26] <br />
- `--perc_het_vars_threshold` :           Proportion of heterozygous variants (%) [15] <br />
- `--scaffold_count_threshold` :          Maximum number of scaffolds per assembly [286] <br />
- `--scaffold_N50_threshold` :            Minimum scaffold N50 [24454] <br />
- `--MASH_hit` :                          Closest MASH hit (top 5) [Streptococcus pneumoniae] <br />

### Output
The output will be written to results/$sample_ID/* which will contain the following results files:

1. `aggregated_stats.long.txt` :        Aggregated statistics row ('long') format.
2. `aggregated_stats.wide.txt` :        Aggregated statistics column ('wide') format.
3. `aggregated_plots.pdf` :             Plots of each metric, with relevant cutoffs plotted. 

