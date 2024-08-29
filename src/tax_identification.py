"""
Module to run Kraken2 and Braken2 on imput fastq files
"""
import os
import subprocess
import pandas as pd
import numpy as np

def run_kraken2(fastq:str, kraken2_db:str, working:str, threads:int=1, fastq2=None):
    """Runs kraken2 on input (paired/unpaired) fastq files using the kraken2_db. Outputs to working"""
    
    #get output and report locations
    output = os.path.join(working, 'kraken2', 'output.txt')
    report = os.path.join(working, 'kraken2', 'report.txt')
    
    #run kraken db on paired reads
    if fastq2 is not None:
        subprocess.run(['kraken2',
                        '--db', kraken2_db,
                        '--threads', str(threads),
                        '--output', output,
                        '--report', report,
                        '--paired',
                        fastq, fastq2], check=True)
    else:
        subprocess.run(['kraken2',
                        '--db', kraken2_db,
                        '--threads', str(threads),
                        '--output', output,
                        '--report', report,
                        fastq], check=True)
        
    return report
        
def run_bracken(kraken_report, 
                kraken2_db:str, 
                working:str, 
                read_length:int=35,
                classification_level='S', 
                read_threshold=10,
                threads:int=1):
    """Performs bayesian restimation for abundance estimation on a kraken2 report"""
    
    ##TODO: error catching if files do not exist
    
    bracken_output = os.path.join(working, 'kraken2', 'output.braken')
    
    print('bracken',
                    '-d', kraken2_db,
                    '-i', kraken_report,
                    '-r', read_length,
                    '-l', classification_level,
                    '-t', read_threshold,
                    '-o', bracken_output)
    
    subprocess.run(['bracken',
                    '-d', kraken2_db,
                    '-i', kraken_report,
                    '-r', read_length,
                    '-l', classification_level,
                    '-t', read_threshold,
                    '-o', bracken_output], check=True)
    
    return bracken_output

    
def get_species_taxid(taxid, ncbi_taxa_db, valid_kingdom):
    lineage = ncbi_taxa_db.get_lineage(taxid)
    if bool(set(lineage) & valid_kingdom):
        taxid2rank_dict = ncbi_taxa_db.get_rank(lineage)
        for lineage_taxid in taxid2rank_dict:
            if taxid2rank_dict[lineage_taxid] == 'species':
                return lineage_taxid
    return None