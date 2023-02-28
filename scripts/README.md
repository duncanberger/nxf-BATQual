# Pipeline and accessory scripts

### asm_stats.py

Calculates standard assembly stats (N50, GC%, assembly size etc), outputting a consistent row format for downstream processing. 

### aggregate.py

Merges results of >=1 run of BATQual and outputs tables and figures for QC and further analyses. 

### run_poppunk.sh

(In development) Runs poppunk on the output of BATQual, implemented as a more efficient means of calculating GPSCs compared to running poppunk, per-sample with each execution of the pipeline. 
- **To do:** Add help message, add better input checks, add parameter to set threads and max memory usage, add final command to delete unnecessary outputs. Add comments. 
