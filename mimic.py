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

from src.tax_identification import run_bracken, run_kraken2

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
        
        kraken = os.path.join(working, 'kraken2')
        os.mkdir(kraken)
        
        magnet = os.path.join(working, 'magnet')
        os.mkdir(magnet)
        
        print('Initialized Working Directory')
    
def run_mimic(args):
    """
    Main pipeline function for the MIMIC
    """
    fastq1 = args.fastq
    fastq2 = args.fastq2
    output = args.output 
    threads = args.threads
    kraken_db = args.db
    
    
    initialize_working(output)
    report_loc = run_kraken2(fastq1, kraken_db, output, threads, fastq2)
    
    ##TODO: Generate braken DB based off read length and kmer-length for the kraken databse (currently pre-generated)
    # report_loc = run_bracken(report_loc, kraken_db, output, threads=threads)
    
    magnet_out = os.path.join(output, 'magnet')
    
    print('python', 'magnet/magnet.py',
                        '-c', report_loc,
                        '-i', fastq1,
                        '-I', fastq2, 
                        '-o', magnet_out,
                        '-t', '4',
                        '-a', '0',
                        '--min-abundance', '0.03')
    if fastq2 is not None:
        subprocess.run(['python', 'magnet/magnet.py',
                        '-c', report_loc,
                        '-i', fastq1,
                        '-I', fastq2, 
                        '-o', magnet_out,
                        '-t', '4',
                        '-a', '0',
                        '--min-abundance', '0.05'], check=True)
    else:
        subprocess.run(['python', 'magnet/magnet.py'
                        '-c', report_loc,
                        '-i', fastq1,
                        '-o', magnet_out,
                        '-t', '4',
                        '-a', '0',
                        '--min-abundance', '0.001'], check=True)
    
    
    
def parse_args():
    parser = argparse.ArgumentParser(description="Universal Taxonomic Classification Verifier.")
    parser.add_argument("-i", "--fastq", type=pathlib.Path, required=True, help="Path to first fastq file")
    parser.add_argument("-I", "--fastq2", type=pathlib.Path, required=False, help="Path to second fastq file for paired-end reads")
    parser.add_argument("-o", "--output", type=pathlib.Path, required=True, help="Path to the output directory.")
    parser.add_argument("--db", type=str, required=True, help='Kraken2 database location')
    parser.add_argument('-t', '--threads', type=int, required=False, default=1, help='Number of threads for multithreading (Default: 1)')


    ## TODO: Decide on other magnet arguments that we can use:
    # parser.add_argument("-m", "--mode", type=str, required=False, choices=['ont', 'illumina'], help="Modes for different sequencing platforms [ont, illumina]. Default:[ont]",  default='ont')
    # parser.add_argument("-t", "--taxid-idx", type=int, required=False, help="The column index (0-based) of the taxids. Default:[0]", default=0)
    # parser.add_argument("-a", "--abundance-idx", type=int, required=False, help="The column index (0-based) of the abundance. Default:[None]")
    # parser.add_argument("--min-abundance", type=float, required=False, help="Minimum abundance (0-1) for pre-filtering, exclude taxa below the threshold.", default=0)
    # parser.add_argument("--min-mapq", type=int, required=False, help="Minimum MAPQ for primary alignments. Default:[20]", default=20)
    # parser.add_argument("--min-covscore", type=float, required=False, help="Minimum Coverage Score for supplementary alignments. Default:[0.7]", default=0.7)
    # parser.add_argument("--kingdom", type=str, help="A comma separated list of taxids of valid kingdoms. Default:[2,4751,2157,10239]", default='2,4751,2157,10239')
    # parser.add_argument("--include-mag", action='store_true', required=False, help="Include metagenomic assemble genomes. Default:[off]")
    # parser.set_defaults(include_mag=False)
    # parser.add_argument("--subspecies", action='store_true', required=False, help="Verify taxonomic classification at subspecies rank. Default:[off]")
    # parser.set_defaults(subspecies=False)
    # parser.add_argument("--accession", action='store_true', required=False, help="Take accession ids as taxids. Does not work with min-abundance. Default:[off]")
    # parser.set_defaults(accession=False)
    
    args = parser.parse_args() 
    
    run_mimic(args) 
    
    

if __name__=='__main__':
    print_info()
    parse_args()
    