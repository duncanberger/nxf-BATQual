# Pipeline and accessory scripts

### asm_stats.py

Calculates standard assembly stats (N50, GC%, assembly size etc), outputting a consistent row format for downstream processing. 

### aggregate.py

Merges results of >=1 run of BATQual and outputs tables and figures for QC and further analyses. 

### run_poppunk.sh

(In development, not implemented) Runs poppunk on the output of BATQual, implemented as a more efficient means of calculating GPSCs compared to running poppunk, per-sample with each execution of the pipeline. 

### cleanup.sh

(In development, not implemented) Simple script designed to remove old work/* files which are no longer needed by BATQual - to save space. 
