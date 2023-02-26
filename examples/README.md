# Testing BATQual

## Downloading test data
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
## Running the pipeline
```
# Run the pipeline using all samples
./nextflow BATQual.nf --input example_fastq.csv --mode fastq --min_k 91 --max_k 133 --outdir example_run
```
## Aggregating results
``` 
# Aggregate all runs
python scripts/aggregate.py --input example_run
```
## Interpreting results

- The pipeline and 'aggregate.py' script will output several files in the 'results/aggregate/' folder. We'll use the 'aggregated_stats.wide.csv' file to flag problematic samples. As a reminder we're expecting:
  -  Two high-quality *S. pnuemoniae* assemblies (pass all checks).
  -  Two low-quality *S. pnuemoniae* assemblies (fail some assembly-quality checks).
  -  One *S. pnuemoniae* assembly contaminated with the same species (Fail species checks and most likely other checks too).
  -  One *S. pnuemoniae* assembly contaminated with a species from the Streptococcus genus (Fail species checks and most likely other checks too).
  -  One *S. pnuemoniae* assembly contaminated with a species from a different genus (Fail species checks and most likely other checks too).
  -  Two Streptococcal assemblies, which are not *S. pnuemoniae* (Fail species checks).
<br/>

### Contiguity $\mathsf{\color{red}{}}$
|Sample|Assembly length (Mb)|Contig N50 (bp)|Contigs (count)|Scaffold N50 (bp)|Scaffolds (count)|
|------|------------------|-------------|------------|---------------|---------------|
|**SAMD00110690**|$\mathsf{\color{red}{4.24}}$|$\mathsf{\color{red}{628}}$|$\mathsf{\color{red}{7492}}$|$\mathsf{\color{red}{634}}$|$\mathsf{\color{red}{7452}}$|
|**SAMEA1023762**|2.05|14803|342|45786|227|
|**SAMEA1024102**|2.17|15818|557|28507|$\mathsf{\color{red}{437}}$|
|**SAMEA1024779**|1.95|45222|194|56344|180|
|**SAMEA1025813**|$\mathsf{\color{red}{5.93}}$|$\mathsf{\color{red}{3314}}$|$\mathsf{\color{red}{5847}}$|10710|$\mathsf{\color{red}{1687}}$|
|**SAMEA1408274**|2.24|41157|221|69004|178|
|**SAMEA2234452**|2.17|54605|152|64127|133|
|**SAMEA2382970**|2.22|66028|255|95094|191|
|**SAMEA3389673**|2.21|74467|323|222119|272|
|**SAMN00761799**|2.07|117890|101|149613|79|
|**SAMN10131018**|$\mathsf{\color{red}{2.28}}$|75456|204|75716|170|
<br/>

- So, we can see that three assemblies have are larger than the default assembly size cutoff (2.26 Mb) and two (SAMD00110690, SAMEA1025813) are far bigger than expected. So these are most likely contaminted with at least one other distinct genome. 
- Both of these assemblies are also very poorly assembled, which also could be a result of contamination preventing proper assembly.
- SAMN10131018 is only marginally over the assembly size limit, so probably something has gone wrong, but it may not be a clear/obvious contaminant. 
<br/>

### Completeness 
|Sample|BUSCO (Complete & single-copy)|BUSCO (Complete & duplicated)|BUSCO (Fragmented)|BUSCO (Missing)|Completeness (CheckM)|Gaps (count)|Gaps (% of assembly)|Gaps sum, bp)
|------|------------------------------|-----------------------------|------------------|---------------|---------------------|------------|--------------------|-------|
|**SAMD00110690**|$\mathsf{\color{red}{29.1}}$|$\mathsf{\color{red}{15.9}}$|$\mathsf{\color{red}{34.1}}$|$\mathsf{\color{red}{20.9}}$|90.1|40|0.03|1082
|**SAMEA1023762**|97.3|0|$\mathsf{\color{red}{1.7}}$|$\mathsf{\color{red}{1}}$|98.44|$\mathsf{\color{red}{115}}$|0.34|6881
|**SAMEA1024102**|96.8|0|$\mathsf{\color{red}{2.7}}$|0.5|98.27|$\mathsf{\color{red}{120}}$|0.26|5561
|**SAMEA1024779**|100|0|0|0|99.58|14|0.03|494
|**SAMEA1025813**|$\mathsf{\color{red}{72.4}}$|$\mathsf{\color{red}{27.4}}$|0.2|0|99.79|$\mathsf{\color{red}{4160}}$|$\mathsf{\color{red}{3.58}}$|$\mathsf{\color{red}{212399}}$
|**SAMEA1408274**|100|0|0|0|99.36|43|0.11|2560
|**SAMEA2234452**|100|0|0|0|99.31|19|0.05|1023
|**SAMEA2382970**|99|0.2|0.2|0.6|$\mathsf{\color{red}{84.56}}$|64|0.17|3824
|**SAMEA3389673**|100|0|0|0|99.62|51|0.37|8106
|**SAMN00761799**|99.8|0|0.2|0|92.74|22|0.05|1097
|**SAMN10131018**|100|0|0|0|93.29|34|0.03|599
<br/>

- We can see that again two assemblies (SAMD00110690, SAMEA1025813) appear to only contain a subset of the expected Benchmarking Universal Single-Copy Orthologs (BUSCO) we'd expect with >15% duplicated again suggesting contaminants in these assemblies. 
- We can also see that two other assemblies (SAMEA1023762, SAMEA1024102) have a slightly higher rate of fragmented BUSCOs than we see in other assemblies and more gaps, suggesting these assemblies are of slightly lower quality than the rest. 
- SAMEA2382970 appears to have the lowest CheckM completeness score but a perfect BUSCO score. The could be a result of variation in the screening gene sets (or screening method) which might be less well suited to this assembly. Which suggests either it's from a different species or there is some issue with the per-base quality which is only flagged by CheckM (meaning it's most likely borderline). 
<br/>

### Species identification 
|Sample|MASH_hit|kraken2_genus|kraken2_species|mlst_species|GC_perc|
|------|--------|-------------|---------------|------------|-------|
|**SAMD00110690**|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.61|
|**SAMEA1023762**|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.63|
|**SAMEA1024102**|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.58|
|**SAMEA1024779**|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.81|
|**SAMEA1025813**|$\mathsf{\color{red}{Bacillus }}$ $\mathsf{\color{red}{ subtilis;}}$ $\mathsf{\color{red}{Alkalihalobacillus }}$ $\mathsf{\color{red}{ gibsonii;}}$ $\mathsf{\color{red}{Bacillus }}$ $\mathsf{\color{red}{ sp.;}}$ $\mathsf{\color{red}{Bacillus }}$ $\mathsf{\color{red}{ subtilis}}$|Streptococcus|Streptococcus pneumoniae|$\mathsf{\color{red}{NA}}$| $\mathsf{\color{red}{42.5}}$|
|**SAMEA1408274**|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.52|
|**SAMEA2234452**|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.61|
|**SAMEA2382970**|$\mathsf{\color{red}{Streptococcus }}$ $\mathsf{\color{red}{ pneumoniae;}}$ $\mathsf{\color{red}{Streptococcus }}$ $\mathsf{\color{red}{ salivarius;}}$ $\mathsf{\color{red}{Streptococcus }}$ $\mathsf{\color{red}{ sp.}}$|Streptococcus|$\mathsf{\color{red}{Streptococcus }}$ $\mathsf{\color{red}{salivarius}}$ |$\mathsf{\color{red}{NA}}$|39.79|
|**SAMEA3389673**|Streptococcus pneumoniae|Streptococcus|Streptococcus pneumoniae|spneumoniae|39.68|
|**SAMN00761799**|$\mathsf{\color{red}{Streptococcus }}$ $\mathsf{\color{red}{ mitis}}$|Streptococcus|$\mathsf{\color{red}{Streptococcus }}$ $\mathsf{\color{red}{ mitis}}$|spneumoniae|39.85|
|**SAMN10131018**|$\mathsf{\color{red}{Streptococcus }}$ $\mathsf{\color{red}{ pseudopneumoniae}}$|Streptococcus|$\mathsf{\color{red}{Streptococcus }}$ $\mathsf{\color{red}{ pseudopneumoniae}}$|spneumoniae|39.91|
<br/>

- SAMEA1025813 was identified primarily as *Bacillus subtilis* by Mash but primarily *S. pnuemoniae* by Kraken2. This would suggest a mixture, probably biased towards *B. subtilis* given the GC% (*B. subtilis* expected GC = 43%) and genome size (expected genome size = 4.21 Mb).
- Despite the oversized assembly, SAMD00110690 appears to be *S. pneumoniae*, both in Mash/Kraken2 hits and GC%. 
- SAMEA2382970 is matches both *S. pneumoniae* and *S. salivarius*. However, there was no evidence of contamination, so most likely it is *S. salivarius*  and this species is just poorly represented in the Mash database, compared to *S. pneumoniae*. 
- SAMN00761799 appears to be *S. mitis*.
- SAMN10131018 appears to be *S. pseudopneumoniae* which explains the slightly larger assembly size. 
<br/>

### Contamination 
|Sample|Heterozygous variants (%)|BUSCO (Complete & duplicated)|Contamination (CheckM)|Strain Heterogenity (CheckM)|GC (%)|
|------|-------------|-------------------------|--------------------|---------------------------|-------|
|**SAMD00110690**|$\mathsf{\color{red}{79.83}}$|$\mathsf{\color{red}{15.9}}$|$\mathsf{\color{red}{117.59}}$|$\mathsf{\color{red}{83.5}}$|39.61|
|**SAMEA1023762**|4.67|0|0.78|0|39.63|
|**SAMEA1024102**|5.46|0|0.38|$\mathsf{\color{red}{50}}$|39.58|
|**SAMEA1024779**|4.82|0|0.83|33.33|39.81|
|**SAMEA1025813**|4.68|$\mathsf{\color{red}{27.4}}$|$\mathsf{\color{red}{85.95}}$|1.63|$\mathsf{\color{red}{42.5}}$|
|**SAMEA1408274**|4.68|0|0.33|0|39.52|
|**SAMEA2234452**|3.89|0|1.33|0|39.61|
|**SAMEA2382970**|$\mathsf{\color{red}{81.60}}$|0.2|$\mathsf{\color{red}{8.24}}$|2.99|39.79|
|**SAMEA3389673**|4.54|0|2.1|0|39.68|
|**SAMN00761799**|0.62|0|2.61|0|39.85|
|**SAMN10131018**|1.96|0|4.15|0|39.91|
<br/>

- SAMD00110690 has a high rate of heterozygous variants suggesting multiple genomes within the assembly and very high contamination and Strain Heterogenity scores consistent with same/different species contamination  and same species contamination, respectively. 
- SAMEA1025813 has a low heterozygous variant rate, but high BUSCO and CheckM contamination scores, which makes sense given that we don't expect *B. subtilis* reads to map to the *S. pneumoniae* assembly.
- SAMEA2382970 does appear to have a very high rate of heterozygous variants but limited evidence of other contamination, suggesting either low level contamination or *S. salivarius* specific nuances (e.g. low mapping rates to the *S. pneumoniae* reference genome resulting in a bias towards heterozygous calls). 
<br/>

### Typing
|Sample|aroE|ddl|gdh|gki|recP|spi|xpt|Serotype (PneumoKITy)|
|------|----|---|---|---|----|---|---|---------------------|
|**SAMD00110690**|588?|863?|$\mathsf{\color{red}{4;11}}$|673?|549?|$\mathsf{\color{red}{4;6}}$|1079?|10A|
|**SAMEA1023762**|10|17|16|40|1|17|1|7F/7A|
|**SAMEA1024102**|1|8|5|4|16|11?|1|14|
|**SAMEA1024779**|15|262|15|34|16|6|14|8|
|**SAMEA1025813**|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|
|**SAMEA1408274**|5|4|6|1|2|6|3|Serogroup_6_(6E)|
|**SAMEA2234452**|5|27|7|4|10|10|1|6A/6B|
|**SAMEA2382970**|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{NA}}$|$\mathsf{\color{red}{Below20-hit-possibleacapsularorganism}}$|
|**SAMEA3389673**|7|8|13|8|6|6|6|35B/35D|
|**SAMN00761799**|~424|142|~477|501|93|393?|~153|$\mathsf{\color{red}{Below70-hit-PoorSequencequality}}$|
|**SAMN10131018**|427|~447|~478|578?|373|~442|~735|$\mathsf{\color{red}{Below20-hit-possibleacapsularorganism}}$|
<br/>

- SAMD00110690 contains multiple MLST alleles, indicative of same species contamination.
- SAMEA1025813 contains no MLST hits, suggesting it's not *S. pneumoniae*.
- SAMEA2382970 contains no MLST hits, suggesting it's not *S. pneumoniae*, but it does have a weak hit to a serotype, indicating it is a different Streptococcal species. 
- SAMN00761799 contains partial MLST hits (?) and novel alleles (~), suggesting it's not *S. pneumoniae*, but it does have a weak hit to a serotype, indicating it is a different Streptococcal species. 
- SAMN10131018 contains partial MLST hits (?) and novel alleles (~), suggesting it's not *S. pneumoniae*, but it does have a weak hit to a serotype, indicating it is a different Streptococcal species. 

### Overall assessment

- **High quality assemblies:** SAMEA1408274, SAMEA2234452
- **Low-quality assemblies:** SAMEA1023762, SAMEA1024102
- **Multiple *S. pneumoniae* genomes:** SAMD00110690
- ***S. salivarius*, possibly also containing *S. pneumoniae*:** SAMEA2382970
- ***S. mitis*:** SAMN10131018
- ***S. pseudopneumoniae*:** SAMN00761799
- ***B. subtilis* and amost certainly some *S. pneumoniae*:** SAMEA1025813
