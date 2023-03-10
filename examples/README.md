# Testing BATQual

Below are a 9 test sequence read sets, covering both high-quality (ideal scenario) inputs and low-quality, contaminated or otherwise problematic inputs. 

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
  -  Two high-quality *S. pneumoniae* assemblies (pass all checks).
  -  Two low-quality *S. pneumoniae* assemblies (fail some assembly-quality checks).
  -  One *S. pneumoniae* assembly contaminated with the same species (Fail species checks and most likely other checks too).
  -  One *S. pneumoniae* assembly contaminated with a species from the Streptococcus genus (Fail species checks and most likely other checks too).
  -  One *S. pneumoniae* assembly contaminated with a species from a different genus (Fail species checks and most likely other checks too).
  -  Two Streptococcal assemblies, which are not *S. pneumoniae* (Fail species checks).
<br/>

### Contiguity 
![image](https://user-images.githubusercontent.com/29282405/224370800-5a2aefe1-c86f-4ca9-9ee6-4f8b7333b827.png)
<br/>

- So, we can see that three assemblies have are larger than the default assembly size cutoff (2.26 Mb) and two (SAMD00110690, SAMEA1025813) are far bigger than expected. So these are most likely contaminted with at least one other distinct genome. 
- Both of these assemblies are also very poorly assembled, which also could be a result of contamination preventing proper assembly.
- SAMN10131018 is only marginally over the assembly size limit, so probably something has gone wrong, but it may not be a clear/obvious contaminant. 
<br/>

### Completeness 
![image](https://user-images.githubusercontent.com/29282405/224370706-214b8806-30dd-4e2e-9fee-2cf35db707eb.png)
<br/>

- We can see that again two assemblies (SAMD00110690, SAMEA1025813) appear to only contain a subset of the expected Benchmarking Universal Single-Copy Orthologs (BUSCO) we'd expect with >15% duplicated again suggesting contaminants in these assemblies. 
- We can also see that two other assemblies (SAMEA1023762, SAMEA1024102) have a slightly higher rate of fragmented BUSCOs than we see in other assemblies and more gaps, suggesting these assemblies are of slightly lower quality than the rest. 
- SAMEA2382970 appears to have the lowest CheckM completeness score but a perfect BUSCO score. The could be a result of variation in the screening gene sets (or screening method) which might be less well suited to this assembly. Which suggests either it's from a different species or there is some issue with the per-base quality which is only flagged by CheckM (meaning it's most likely borderline). 
<br/>

### Species identification 
![image](https://user-images.githubusercontent.com/29282405/224370486-88826f64-7412-4e8c-b9e4-7562cc4e0842.png)
<br/>

- SAMEA1025813 was identified primarily as *Bacillus subtilis* by Mash but primarily *S. pneumoniae* by Kraken2. This would suggest a mixture, probably biased towards *B. subtilis* given the GC% (*B. subtilis* expected GC = 43%) and genome size (expected genome size = 4.21 Mb).
- Despite the oversized assembly, SAMD00110690 appears to be *S. pneumoniae*, both in Mash/Kraken2 hits and GC%. 
- SAMEA2382970 is matches both *S. pneumoniae* and *S. salivarius*. However, there was no evidence of contamination, so most likely it is *S. salivarius*  and this species is just poorly represented in the Mash database, compared to *S. pneumoniae*. 
- SAMN00761799 appears to be *S. mitis*.
- SAMN10131018 appears to be *S. pseudopneumoniae* which explains the slightly larger assembly size. 
<br/>

### Contamination 
![image](https://user-images.githubusercontent.com/29282405/224370221-7793be1d-5f18-4d1c-8d4f-b082783bd921.png)
<br/>

- SAMD00110690 has a high rate of heterozygous variants suggesting multiple genomes within the assembly and very high contamination and Strain Heterogenity scores consistent with same/different species contamination  and same species contamination, respectively. 
- SAMEA1025813 has a low heterozygous variant rate, but high BUSCO and CheckM contamination scores, which makes sense given that we don't expect *B. subtilis* reads to map to the *S. pneumoniae* assembly.
- SAMEA2382970 does appear to have a very high rate of heterozygous variants but limited evidence of other contamination, suggesting either low level contamination or *S. salivarius* specific nuances (e.g. low mapping rates to the *S. pneumoniae* reference genome resulting in a bias towards heterozygous calls). 
<br/>

### Typing
![image](https://user-images.githubusercontent.com/29282405/224370078-2ffb37ee-1632-4a5b-bfaf-74576da0c222.png)

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
