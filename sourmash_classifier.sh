#!/bin/bash

# Make a directory for the kmers
mkdir -p ~/kmers
cd ~/kmers || exit # Change to kmers directory, exit if fails

# Define the location for the sourmash database and taxanomy spreadsheet
sourmash_db=~/home/sourmash_database
taxa_spreadsheet=~/home/sourmash_script

# Index the taxanomy spreadsheet otherwise provide the indexed spreadsheet for taxa spreadsheet
sourmash tax prepare -t "$taxa_spreadsheet" -o "$(basename $taxa_spreadsheet).taxonomy.sqldb" -F sql

# Define path for the input reads
input_reads=~/samples/path/  

# Ensure the directory exists
if [ ! -d "$input_reads" ]; then
    echo "Input reads directory does not exist: $input_reads"
    exit 1
fi

# Build signature files
for file in "$input_reads"/*; do
    if [ -f "$file" ]; then
        output_file=$(basename "$file" | cut -f1 -d"_")
        sourmash sketch dna -p k=31,abund "$file" -o "${output_file}.sig.gz" --name "$output_file"
    fi
done

# Find matching genomes with Sourmash gather
for sig in *.sig.gz; do
    sourmash gather "$sig" $(basename "$sourmash_db" .zip) --save-matches matches.zip
    
    # Rerun gather, save the results to a CSV
    sourmash gather "$sig" matches.zip -o "${sig%.sig.gz}.x.gtdb.csv"
done

# Use tax metagenome to classify the metagenome
echo "Classifying the metagenome using tax metagenome..."
for csv in *.x.gtdb.csv; do
    sourmash tax metagenome -g "$csv" -t taxonomy.sqldb -F human -r order -o "${csv%.csv}_classified.csv"
done

# End