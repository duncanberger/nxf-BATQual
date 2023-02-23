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

# Download contaminated (same species) Streptococcus pneumoniae FASTQ files (SAMEA1025813) 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/007/ERR1202187/ERR1202187_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/007/ERR1202187/ERR1202187_2.fastq.gz

# Download Streptococcus pneumoniae FASTQ files that yield low quality assemblies (SAMEA1023762)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047984/ERR047984_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047984/ERR047984_2.fastq.gz

# Download Streptococcus pneumoniae FASTQ files that yield low quality assemblies (SAMEA1024779)
 wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050026/ERR050026_1.fastq.gz
 wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050026/ERR050026_2.fastq.gz

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

### Aggregating results

### Checking results

