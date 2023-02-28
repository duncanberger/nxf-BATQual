#!/bin/bash

# Define default values
sample_list=""
baseDir=""

# Parse command-line parameters
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        --sample_list)
            sample_list="$2"
            shift # past argument
            shift # past value
            ;;
        --baseDir)
            baseDir="$2"
            shift # past argument
            shift # past value
            ;;
        *)    # unknown option
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if required parameters are present
if [ -z "$sample_list" ] || [ -z "$baseDir" ]
then
    echo "Missing required parameter(s). Usage: $0 --sample_list SAMPLE_LIST --baseDir BASE_DIR"
    exit 1
fi

# Execute code for each sample in the list
while read -r sample_id
do
    echo $sample_id $baseDir/$sample_id/$sample_id.velvet_contigs.fa | tr ' ' '\t'
done < "$sample_list" > run_GPSC.list

# Run Poppunk on the combined input (giving a list of all relevant accessions)
poppunk_assign --db $baseDir/DB/GPS_v6 --distances $baseDir/DB/GPS_v6/GPS_v6.dists --query GPSC_input.txt --output out_GPSC --external-clustering $baseDir/DB/GPS_v6_external_clusters.csv --threads 5

# Split the result into individual files and place them back in each samples directory
while read -r sample_id
do
grep $sample_id out_GPSC_external_clusters.csv | awk -F, '{print $1,"GPSC",$2,""}' OFS=',' > $baseDir/$sample_id/$sample_id.GPSC_results.txt
done < "$sample_list"
