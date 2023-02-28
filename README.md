# BATQual [Pre-release]
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.4-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)

## Table of contents 
* [Introduction](#Introduction)
* [Pipeline summary](#pipeline_summary)
* [Installation](#install)
* [Running BATQual](#run)
* [Aggregation and quality control](#aggregate)
* [Run PopPUNK (optional)](#poppunk)
* [Example datasets](#examples)
* [Tips: Memory use and efficiency](#mem)
* [Tips: Adjusting config files](#tips_adj)
* [Components](#components)
* [Citation](#cite)


## Introduction <a name="Introduction"></a>
**BATQual** is a Nextflow pipeline for performing assembling, annotation and typing bacterial genomes. It's designed to work on either Illumina short-reads or assembled genomes. It was primarily developed to assemble and type _Streptococcus pneumoniae_ genomes and so there is a specific mode '--pneumo' which will enable serotyping and cluster assignment tools specific to this species, if it's required. 

**[Note: BATQual is currently being developed and has only been run internally]** <br />

**[Note: Addition of executors to the pipeline is in progress - most likely will default to Slurm]**

**[Note: Conda environment is functional but needs to be simplified]**


## Pipeline summary <a name="pipeline_summary"></a>

### Schematic overview 
![rect13856](https://user-images.githubusercontent.com/29282405/219794605-086076a0-2e6d-471b-8e48-501d00f853b2.png)

### Pipeline description
#### FASTQ processing
When processing Illumina sequencing reads ('--mode fastq'), the pipeline will perform the following steps:

1. Trim adaptors from reads ([`fastp`](https://github.com/OpenGene/fastp)).
2. Assess read quality using ([`FASTQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)). 
3. Assemble one genome per-isolate using ([`VelvetOptimiser`](https://github.com/tseemann/VelvetOptimiser)).
4. Estimate heterozygosity with the read data to check for contamination using ([`minimap`](https://github.com/lh3/minimap2)) and ([`BCFtools`](https://samtools.github.io/bcftools/bcftools.html)).
5. Estimate read contamination using([`Kraken2`](https://github.com/DerrickWood/kraken2))
6. Calculate standard assembly statistics using a custom script ([`asm_stats.py`](https://github.com/duncanberger/nxf-BATQual/blob/main/scripts/asm_stats.py))
7. Estimate assembly completeness using ([`BUSCO`](https://busco.ezlab.org/))
8. Check for assembly completeness and contamination using ([`CheckM`](https://github.com/Ecogenomics/CheckM))
9. Evaluate assembly using ([`QUAST`](https://quast.sourceforge.net/)) 
10. Estimate most likely species based on the closet ([`MASH`](https://github.com/marbl/Mash)) hits in the RefSeq database
11. Perform read-based pneumococcal serotyping using ([`seroBA`](https://github.com/sanger-pathogens/seroba))
12. Perform assembly-based pneumococcal serotyping using ([`pneumoKITy`](https://github.com/sanger-pathogens/seroba))
13. Assign Global Pneumococcal Sequence Clusters ([`GPSCs`](https://www.pneumogen.net/gps/)) using popunk ([`PopPunk`](https://poppunk.net/)) 
14. Perform multilocus sequence typing using ([`mlst`](https://github.com/tseemann/mlst))
15. Annotate each genome using ([`Prokka`](https://github.com/tseemann/prokka)) 

#### FASTA processing
When processing genome assemblies ('--mode fasta'), the pipeline will omit read-specific steps:
1. Calculate standard assembly statistics using a custom script ([`asm_stats.py`](https://github.com/duncanberger/nxf-BATQual/blob/main/scripts/asm_stats.py))
2. Estimate assembly completeness using ([`BUSCO`](https://busco.ezlab.org/))
3. Check for assembly completeness and contamination using ([`CheckM`](https://github.com/Ecogenomics/CheckM))
4. Evaluate assembly using ([`QUAST`](https://quast.sourceforge.net/)) 
5. Estimate most likely species based on the closet ([`MASH`](https://github.com/marbl/Mash)) hits in the RefSeq database
6. Perform assembly-based pneumococcal serotyping using ([`pneumoKITy`](https://github.com/sanger-pathogens/seroba))
7. Assign Global Pneumococcal Sequence Clusters ([`GPSCs`](https://www.pneumogen.net/gps/)) using popunk ([`PopPunk`](https://poppunk.net/)) 
8. Perform multilocus sequence typing using ([`mlst_check`](https://github.com/sanger-pathogens/mlst_check))
9. Annotate each genome using ([`Prokka`](https://github.com/tseemann/prokka)) 

## Installation <a name="install"></a>
### Software
You will need to install [`Conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and/or [`Mamba`](https://github.com/conda-forge/miniforge#mambaforge). <br />

You will need to install [`pneumoKITy`](https://github.com/CarmenSheppard/PneumoKITy), instructions can be found on the PneumoKITy github. 

### Installing BATQual
```
# Clone this repo
git clone https://github.com/duncanberger/nxf-BATQual.git

# Set up conda dependencies 
conda env create -f environment.yml 

# Activate environment
conda activate nxf-BATQual
```
### Databases
This pipeline uses a number of databases, before running the pipeline you should check that version included below are the most relevant/up to date for your analyses.  The can be downloaded and made ready for analysis as follows:

```
# If not already present make a folder for your databases in the nxf-BATQual folder
mkdir -p DB/
```

#### RefSeq Mash database
```
# Download a preformatted MASH database of RefSeq (release 70, if you want a more up to date version you'll need to build your own)
wget https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh

# Move it to the database folder
mv refseq.genomes.k21s1000.msh DB/

# Download a file containing database IDs and species'
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt

# Move it to the database folder
mv assembly_summary_genbank.txt DB/

```
#### seroBA
```
# seroBA can be found here: https://github.com/sanger-pathogens/seroba.git but all we want is the database folder. You can either clone the whole repo, and grab the database folder (and delete the rest):
git clone https://github.com/sanger-pathogens/seroba.git

# Or you can download that folder specifically:
git init
git remote remove origin
git remote add origin -f https://github.com/sanger-pathogens/seroba.git
echo "database" > .git/info/sparse-checkout
git config core.sparseCheckout true
git pull origin master

# In either case just move the database folder into the DB/ folder:
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

# Index the FASTA
cd DB/
samtools index GCA_000026665.1_ASM2666v1_genomic.fna

# Then add the assembly name to 'main.config' file at the 'ref_assembly=' row
```

#### Benchmarking Universal Single-Copy Orthologs (BUSCO)
```
# By default BUSCO will download the relevant databases when trying to run, however, when running in nextflow it'll try to do that for every sample. A better solution is to download the relevant database prior to running the pipeline and specifying the relevant path. 
wget https://busco-data.ezlab.org/v5/data/lineages/lactobacillales_odb10.2020-03-06.tar.gz
mv lactobacillales_odb10.2020-03-06.tar.gz DB/
tar xvf lactobacillales_odb10.2020-03-06.tar.gz
rm lactobacillales_odb10.2020-03-06.tar.gz
```
#### Kraken 2 database
```
# Download the Kraken DB (capped at 8 Gb):
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20221209.tar.gz

# Move to database directory and unpack
mv k2_standard_08gb_20221209.tar.gz DB/
cd DB/
tar xvf k2_standard_08gb_20221209.tar.gz
```

## Running BATQual <a name="run"></a>

You can run the pipeline as follows:
```
nextflow run BATQual.nf --input fastq_files.csv --mode fastq
```
or 
```
nextflow run BATQual.nf --input fasta_files.csv --mode fasta
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

1. `$sample_ID.velvet_contigs.fa` : VelvetOptimiser genome assembly
2. `$sample_ID.prokka.fna` : VelvetOptimiser genome assembly renamed by Prokka to match GFF file
3. `$sample_ID.prokka.gff` : VelvetOptimiser genome annotation, in Prokka format
4. `$sample_ID.prokka_results.txt` : Counts of predicted genes
5. `$sample_ID.hetperc_results.txt` : Estimated heterozygosity of sequencing reads
6. `$sample_ID.assembly_results.txt` : Assembly statistics
7. `$sample_ID.BUSCO_results.txt` : seroBA serotyping results 
8. `$sample_ID.checkM_results.txt` : seroBA serotyping results
9. `$sample_ID.mash_results.txt` : seroBA serotyping results
10. `$sample_ID.seroBA_results.txt` : seroBA serotyping results
11. `$sample_ID.pneumoKITy_results.txt` : seroBA serotyping results
12. `$sample_ID.GPSC_results.txt` : seroBA serotyping results 
13. `$sample_ID.kraken2_results.txt` : Kraken2 read screening results 
14. `$sample_ID.kraken2.out` : Raw Kraken2 output
14. `$sample_ID.R*_fastqc.zipt` : Raw FASTQC results 
14. `$sample_ID.quast` : QUAST output directory


Each file is formatted into columns as follows: 
1. `sample_ID`: Isolate ID
2. `metric`: Metric used for evaluation
3. `results`: Metric results
4. `status`: QC status, unpopulated except for pneumoKITy QC, to be populated at the next stage. 

## Run PopPUNK (optional) <a name="poppunk"></a>

PopPunk can be run as part of BATQual but as it runs individually (one sample at a time), it nearly doubles the runtime. You could run it as a final step (after all genomes have been assembled) to sidestep this problem, but you may want to run BATQual multiple times across different sample sets and so it's much simpler to split off PopPUNK GPSC assignment into it's own script.
```
# That script can be found and run in the following way:
scripts/run_poppunk_FASTQ.sh --sample_list SAMPLE_LIST --baseDir BASE_DIR --threads THREADS
```
Where 'SAMPLE_LIST' is a single column list of samples and baseDir is the folder where BATQual was run (the location where the output folder was created.

If you're just working with FASTA files then use this script:
```
scripts/run_poppunk_FASTA.sh --input_list INPUT_LIST --baseDir BASE_DIR --threads THREADS
```
Where 'INPUT_LIST' is the same *.csv file you provided to BATQual (col 1:sample_id, col 2: FASTA location) with your FASTA file locations. baseDir is the folder where BATQual was run (the location where the output folder was created.

Once either of these scripts have run, it will output your GPSC assignment results in each sample directory which can then be further processed with 'aggregate.py'.


## Aggregation and quality control <a name="aggregate"></a>
To aggregate the results across multiple isolates and runs, I have written an accessory script to produce merged output tables and plot the results. Again 

You can run the script as follows (where results is the name of the folder specified by the '--output' parameter:
```
python scripts/aggregate.py --input results
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
4. `multiqc_report.html` :              MultQC reports (FASTQC, QUAST and Kraken2 summaries)

## Example datasets <a name="examples"></a>
A full guide to running through the pipeline with exemplar read datasets for a variety of possible scenarios can be found [`here`](https://github.com/duncanberger/nxf-BATQual/tree/main/examples).

## Tips: Runtimes and efficiency <a name="mem"></a>
- You can increase the rate of sample processing in a few ways. The most time consuming process is the Velvet assembly stage this can be sped up by decreasing the range of kmer sizes ('--min_k' and '--max_k' parameters) which VelvetOptimiser uses. Higher read-depth inputs tend to have optimum kmer lengths ~127 bp, while lower read-depth/more error prone datasets are more suited to shorter kmers (70-95 bp). 
- Disabling the more time consuming processes can decrease run times. These include: the het_call, CheckM and Kraken2 stages. The downside being that you will lose these quality checks, but the option exists to turn these steps off. 
- By default GPSCs are not assigned but the process exists and can be enabled. Poppunk is intended to run on multiple samples, but by default is run per-sample in this pipeline (for scenarios where you might only want to assemble a few samples). If you want to run Poppunk GPSC assignment at scale, I'd recommend running 'scripts/run_GPSC.sh' script. This will deposit output files, per-sample, in the matching directories, which can then be processed by the 'scripts/aggregate.py' script.

## Tips: Adjusting config files <a name="tips_adj"></a>
- Parameters for the execution of BATQual can be found in the 'nextflow.config' file. These can be also be set on the command line during execution. 

## Components <a name="components"></a>

Key tools used by BATQual:
- bcftools=1.8
- busco=5.4.5
- biopython=1.77
- checkm-genome=1.2.2
- fastp=0.23.2
- fastqc=0.11.9
- kraken2=2.1.2
- minimap2=2.24
- mlst=2.23.0
- multiqc=1.14
- nextflow=22.10.6
- numpy=1.24.2
- pandas=1.5.3
- pneumoKITy=1.0
- poppunk=2.6.0
- prokka=1.14.6
- python=3.9.15
- quast=5.2.0
- samtools=1.15.1
- scipy=1.10.0
- seaborn=0.12.2
- seroba=1.0.2
- VelvetOptimiser=2.2.6

## Citation <a name="cite"></a>
