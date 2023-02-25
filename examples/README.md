# Testing BATQual

### Downloading test data
```
# Download high-quality Streptococcus pneumoniae FASTQ files (SAMEA1408274)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/ERR107504/ERR107504_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/ERR107504/ERR107504_2.fastq.gz
                                                                                                                                                             
# Download high-quality Streptococcus pneumoniae FASTQ files (SAMEA2234452)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR449/ERR449780/ERR449780_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR449/ERR449780/ERR449780_2.fastq.gz

# Download contaminated (different species) Streptococcus pneumoniae FASTQ files (SAMEA1025813) 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR066/ERR066215/ERR066215_1.fastq.gz 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR066/ERR066215/ERR066215_2.fastq.gz

# Download contaminated (same species) Streptococcus pneumoniae FASTQ files (SAMD00110690)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR121/DRR121429/DRR121429_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR121/DRR121429/DRR121429_2.fastq.gz

# Download Streptococcus pneumoniae FASTQ files that yield low quality assemblies (SAMEA1023762)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047984/ERR047984_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047984/ERR047984_2.fastq.gz

# Download Streptococcus pneumoniae FASTQ files that yield low quality assemblies (SAMEA1024102)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047972/ERR047972_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047972/ERR047972_2.fastq.gz

# Download FASTQ files of a species falsely identified as Streptococcus pneumoniae (SAMEA2382970)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR563/ERR563704/ERR563704_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR563/ERR563704/ERR563704_2.fastq.gz

# Download Streptococcus pseudopneumoniae FASTQ files (SAMN10131018)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/009/SRR7904309/SRR7904309_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/009/SRR7904309/SRR7904309_2.fastq.gz

# Download Streptococcus mitis FASTQ files (SAMN00761799)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR387/SRR387657/SRR387657_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR387/SRR387657/SRR387657_2.fastq.gz

```
### Running the pipeline
```
# Run the pipeline using all samples
./nextflow BATQual.nf --input example_fastq.csv --mode fastq --min_k 91 --max_k 133 --outdir example_run
```
### Aggregating results
``` 
# Aggregate all runs
python scripts/aggregate.py --input example_run
```
### Interpreting results
#### Contiguity
sample|Assembly length (Mb)|Contig N50 (bp)|Contigs (count)|Scaffold N50 (bp)|Scaffolds (count)|
|------|------------------|-------------|------------|---------------|---------------|
|SAMD00110690|$\mathsf{\color{red}{4.24}}$|$\mathsf{\color{red}{628}}$|$\mathsf{\color{red}{7492}}$|$\mathsf{\color{red}{634}}$|$\mathsf{\color{red}{7452}}$|
|SAMEA1023762|2.05|14803|342|45786|227|
|SAMEA1024102|2.17|15818|557|28507|437|
|SAMEA1024779|1.95|45222|194|56344|180|
|SAMEA1025813|5.93|3314|5847|10710|1687|
|SAMEA1408274|2.24|41157|221|69004|178|
|SAMEA2234452|2.17|54605|152|64127|133|
|SAMEA2382970|2.22|66028|255|95094|191|
|SAMEA3389673|2.21|74467|323|222119|272|
|SAMN00761799|2.07|117890|101|149613|79|
|SAMN10131018|2.28|75456|204|75716|170|

#### Completeness

|sample|BUSCO (Complete & single-copy)|BUSCO (Complete & duplicated)|BUSCO (Fragmented)|BUSCO (Missing)|Completeness (CheckM)|Gaps (count)|Gaps (% of assembly)|Gaps sum, bp)
|------|------------------------------|-----------------------------|------------------|---------------|---------------------|------------|--------------------|-------|
SAMD00110690|$\mathsf{\color{red}{29.1}}$|$\mathsf{\color{red}{15.9}}$|$\mathsf{\color{red}{34.1}}$|$\mathsf{\color{red}{20.9}}$|90.1|40|0.03|1082
SAMEA1023762|97.3|0|1.7|1|98.44|115|0.34|6881
SAMEA1024102|96.8|0|2.7|0.5|98.27|120|0.26|5561
SAMEA1024779|100|0|0|0|99.58|14|0.03|494
SAMEA1025813|72.4|27.4|0.2|0|99.79|4160|3.58|212399
SAMEA1408274|100|0|0|0|99.36|43|0.11|2560
SAMEA2234452|100|0|0|0|99.31|19|0.05|1023
SAMEA2382970|99|0.2|0.2|0.6|84.56|64|0.17|3824
SAMEA3389673|100|0|0|0|99.62|51|0.37|8106
SAMN00761799|99.8|0|0.2|0|92.74|22|0.05|1097
SAMN10131018|100|0|0|0|93.29|34|0.03|599

#### Species identification

|sample|MASH_hit|kraken2_genus|kraken2_species|mlst_species|GC_perc|
|------|--------|-------------|---------------|------------|-------|
|SAMD00110690|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.61|
|SAMEA1023762|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.63|
|SAMEA1024102|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.58|
|SAMEA1024779|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.81|
|SAMEA1025813|Bacillus subtilis; Alkalihalobacillus gibsonii; Bacillus sp.; Bacillus subtilis|Streptococcus|Streptococcus pneumoniae|NA|42.5|
|SAMEA1408274|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.52|
|SAMEA2234452|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.61|
|SAMEA2382970|Streptococcus pneumoniae; Streptococcus salivarius; Streptococcus sp.|Streptococcus|Streptococcus salivarius|-|39.79|
|SAMEA3389673|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.68|
|SAMN00761799|Streptococcus mitis|Streptococcus|Streptococcus mitis|spneumoniae|39.85|
|SAMN10131018|Streptococcus pseudopneumoniae|Streptococcus|Streptococcus pseudopneumoniae|spneumoniae|39.91|

#### Contamination $\mathsf{\color{red}{}}$
sample|Heterozygous variants (%)|BUSCO (Complete & duplicated)|Contamination (CheckM)|Strain Heterogenity (CheckM)|GC (%)|
------|-------------|-------------------------|--------------------|---------------------------|-------|
SAMD00110690|$\mathsf{\color{red}{79.83}}$|$\mathsf{\color{red}{15.9}}$|$\mathsf{\color{red}{117.59}}$|$\mathsf{\color{red}{83.5}}$|39.61|
SAMEA1023762|4.67|0|0.78|0|39.63|
SAMEA1024102|5.46|0|0.38|50|39.58|
SAMEA1024779|4.82|0|0.83|33.33|39.81|
SAMEA1025813|4.68|$\mathsf{\color{red}{27.4}}$|$\mathsf{\color{red}{85.95}}$|1.63|$\mathsf{\color{red}{42.5}}$|
SAMEA1408274|4.68|0|0.33|0|39.52|
SAMEA2234452|3.89|0|1.33|0|39.61|
SAMEA2382970|$\mathsf{\color{red}{81.60}}$|0.2|$\mathsf{\color{red}{8.24}}$|2.99|39.79|
SAMEA3389673|4.54|0|2.1|0|39.68|
SAMN00761799|0.62|0|2.61|0|39.85|
SAMN10131018|1.96|0|4.15|0|39.91|

#### Typing
|sample|aroE|ddl|gdh|gki|recP|spi|xpt|Serotype (PneumoKITy)|
|------|----------------|---------------|---------------|---------------|----------------|---------------|---------------|-------------------|
|SAMD00110690|588?|863?|$\mathsf{\color{red}{4;11}}$|673?|549?|$\mathsf{\color{red}{4;6}}$|1079?|10A|
|SAMEA1023762|10|17|16|40|1|17|1|7F/7A|
|SAMEA1024102|1|8|5|4|16|11?|1|14|
|SAMEA1024779|15|262|15|34|16|6|14|8|
|SAMEA1025813|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|
|SAMEA1408274|5|4|6|1|2|6|3|Serogroup_6_(6E)|
|SAMEA2234452|5|27|7|4|10|10|1|6A/6B|
|SAMEA2382970|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{Below20-hit-possibleacapsularorganism}}$|
|SAMEA3389673|7|8|13|8|6|6|6|35B/35D|
|SAMN00761799|~424|142|~477|501|93|393?|~153|$\mathsf{\color{red}{Below70-hit-PoorSequencequality}}$|
|SAMN10131018|427|~447|~478|578?|373|~442|~735|$\mathsf{\color{red}{Below20-hit-possibleacapsularorganism}}$|
