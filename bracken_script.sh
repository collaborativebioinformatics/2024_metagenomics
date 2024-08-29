#!/bin/bash

# This script is for running bracken on the kraken output

# Define the location of the kraken database
kraken_db_loc=~/home/kraken2_database

# Define the directory where kraken2 outputs are stored from step 1
kraken_output_dir=~/home/kraken_output

# Define the directory where bracken outputs will be stored
bracken_output_dir=~/home/bracken_output

# Listing files ending with .fastq.gz, cutting the first 3 fields using underscore as a delimiter,
# and sorting uniquely to handle paired reads or sample prefixes
for i in `ls samples/*fastq.gz | cut -f1-3 -d"_" | sort -u`
do
    # Extract the sample name from the variable $i by splitting on the slash and taking the second field
    output_file=`echo $i | cut -f2 -d"/"`

    # Run Bracken using the Kraken2 database and the corresponding Kraken2 report for the sample
    # Input is the Kraken2 report file and output is the Bracken results file
    bracken -d $kraken_db_loc \
            -i ${kraken_output_dir}/${output_file}.krakenreport \
            -o ${bracken_output_dir}/${output_file}.bracken
done

# End

