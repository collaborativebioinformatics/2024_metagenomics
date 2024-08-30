''' 
MIMIC Software

@Authors:
Ryan Doughty, Iva Kotaskova, Kasambula Arthur Shem, Shwetha Kumar, Mike Nute, Todd Treangen, Mike Nute

@Version 0.1


Mimic creates simulated metagenomes based off of real existing metagenomes.
'''

import argparse
import os
import pathlib
import subprocess

from src.tax_identification import run_bracken, run_kraken2, run_lemur, parse_magnet_output
from src.sim import generate_species_file_info, prep_sim_lemur, run_read_analysis, run_sim

__author__ = "Ryan Doughty, Iva Kotaskova, Kasambula Arthur Shem, Shwetha Kumar, Mike Nute, Todd Treangen"
__contact__ = "rdd4@rice.edu"

__license__ = "MIT"
__version__ = "0.1"
__email__ = "rdd4@rice.edu"
__status__ = "Development"

def print_info():
    """
    Prints tool information
    """
    print('----------------------------------------')
    print(f'MIMIC Software Version {__version__}')
    print('----------------------------------------')
    print(f'Contact: {__contact__}')
    print(f'Authors: {__author__}\n')
    print(f'Status: {__status__}\n\n')
    
def initialize_working(working:str):
    """Initializes working directory, creates necessary sub directories"""

    if os.path.exists(working):
        raise SystemExit('Initializing Working Directory Failed, working directory already exists')
    else:
        os.mkdir(working)
        
        lemur = os.path.join(working, 'lemur')
        os.mkdir(lemur)
        
        magnet = os.path.join(working, 'magnet')
        os.mkdir(magnet)
        
        nanosim = os.path.join(working, 'nanosim')
        os.mkdir(nanosim)
        
        print('Initialized Working Directory\n')
    
def run_mimic(args):
    """
    Main pipeline function for the MIMIC
    """
    fastq1 = args.fastq
    fastq2 = args.fastq2
    output = args.output 
    threads = args.threads
    lemur_db = args.db
    num_reads = args.reads
    
    
    initialize_working(output)
    
    ## run lemur
    lemur_out = os.path.join(output, 'lemur')
    report_loc = run_lemur(fastq1, lemur_db, lemur_out, threads=threads)
    
    magnet_out = os.path.join(output, 'magnet')
    
    if fastq2 is not None:
        subprocess.run(['python', 'magnet/magnet.py',
                        '-c', report_loc,
                        '-i', fastq1,
                        '-I', fastq2, 
                        '-o', magnet_out,
                        '--min-abundance', '0.1',
                        '--threads', str(threads)], check=True)
    else:
        subprocess.run(['python', 'magnet/magnet.py',
                        '-c', report_loc,
                        '-i', fastq1,
                        '-o', magnet_out,
                        '-a', '12',
                        '--min-abundance', '0.01',
                        '--threads', str(threads)], check=True)
    
    magnet_report = os.path.join(magnet_out, 'cluster_representative.csv')
    if not os.path.exists(magnet_report):
        raise SystemExit('Magnet failed')
    
    
    ## get necessary files for the nanosim input and run nanosim
    nanosim_loc = os.path.join(output, 'nanosim')
    genome_list = prep_sim_lemur(magnet_report, report_loc, nanosim_loc, num_reads)
    species_info = generate_species_file_info(genome_list, nanosim_loc)
    
    genome_list_loc = os.path.join(nanosim_loc, 'genome_list1.tsv')
    run_read_analysis(fastq1, genome_list_loc, nanosim_loc, threads=threads) ## nanosim step 1

    genome_list_loc = os.path.join(nanosim_loc, 'genome_list2.tsv')
    abundances = os.path.join(nanosim_loc, 'abundances.tsv')
    species_loc = os.path.join(nanosim_loc, 'species_info.tsv')
   
    run_sim(genome_list_loc, abundances, species_loc, nanosim_loc, threads=threads) ## nanosim step 2

    
def parse_args():
    parser = argparse.ArgumentParser(description="Universal Taxonomic Classification Verifier.")
    parser.add_argument("-i", "--fastq", type=pathlib.Path, required=True, help="Path to first fastq file")
    parser.add_argument("-I", "--fastq2", type=pathlib.Path, required=False, help="Path to second fastq file for paired-end reads")
    parser.add_argument("-o", "--output", type=pathlib.Path, required=True, help="Path to the output directory.")
    parser.add_argument("--db", type=str, required=True, help='Kraken2 database location')
    parser.add_argument('-t', '--threads', type=int, required=False, default=1, help='Number of threads for multithreading (Default: 1)')
    parser.add_argument('-r', '--reads', type=int, required=True, default=100, help='Number of simulated reads to generate')
    
    args = parser.parse_args() 
    
    run_mimic(args) 
    
    

if __name__=='__main__':
    print_info()
    parse_args()
    