# Pipeline and accessory scripts

### asm_stats.py

Calculates standard assembly stats (N50, GC%, assembly size etc), outputting a consistent row format for downstream processing. 

### aggregate.py

Merges results of >=1 run of BATQual and outputs tables and figures for QC and further analyses. 

### run_poppunk_FASTQ.sh
Runs poppunk on the output of BATQual, implemented as a more efficient means of calculating GPSCs compared to running poppunk, per-sample with each execution of the pipeline. 
- **To do:** Add help message.

### run_poppunk_FASTA.sh
Runs poppunk on the output of BATQual, implemented as a more efficient means of calculating GPSCs compared to running poppunk, per-sample with each execution of the pipeline. 
- **To do:** Add help message.
