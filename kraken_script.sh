#!/bin/bash
# This is a script for processing paired end reads using kraken2
# It will loop over all the samples

# Define the location of the kraken2 database
kraken_db_loc=~/kraken2_database

# Define the output directory for the kraken2 results
output_dir=~/kraken_output

# Create the output directory if it doesnot exist
mkdir -p $output_dir

# Loop through all the pairs of fastq.gz files in the sample directory
# This assumes that all samples are in the samples directory and are compressed fastq files with the extension fastq.gz
# It also assumes that they follow a particular naming style

for R1 in samples/*_R1_*.fastq.gz; do
    # Identify the corresponding R2 file by replacing R1 with R2 in the filename
    R2=${R1/_R1_/_R2_}

    # Extract the base name of the file without the directory path and extensions
    output_file=$(basename "$R1" | cut -f1 -d"_")

    # Run Kraken2 on the paired files
    kraken2 --db $kraken_db_loc --paired $R1 $R2 --report $output_dir/${output_file}.krakenreport > $output_dir/${output_file}.kraken2

done
# End of the for loop.
