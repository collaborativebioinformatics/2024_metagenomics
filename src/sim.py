import os
import pandas as pd
import numpy as np
import subprocess

def run_read_analysis(fastq:str, genome_list:str, out_loc:str, threads:int=1):
    
    subprocess.run(['read_analysis.py',
                    'metagenome',
                    '-i', fastq,
                    '-gl', genome_list,
                    '-o', out_loc + '/training/training',
                    '-t', str(threads)], check=True)
    
def run_sim(genome_list:str, abundance_list:str, species_list:str, out_loc:str, perfect:bool=False, threads:int=1):
    
    if perfect:
        subprocess.run(['simulator.py', 'metagenome',
                        '-gl', genome_list,
                        '-a', abundance_list,
                        '-dl', species_list,
                        '-c', out_loc + '/training/training',
                        '-o', out_loc + '/simulated',
                        '--perfect',
                        '-t', str(threads)], check=True)
    else:
        subprocess.run(['simulator.py', 'metagenome',
                        '-gl', genome_list,
                        '-a', abundance_list,
                        '-dl', species_list,
                        '-c', out_loc + '/training/training',
                        '-o', out_loc + '/simulated',
                        '-t', str(threads)], check=True)

def get_lemur_abundance(taxid, lemur_data):
    try:
        out = lemur_data[lemur_data['Target_ID'] == taxid]
        if not out.empty:
            return out['F'].values[0]
        else:
            return 0
    except KeyError:
        return 0

def prep_sim_lemur(metadata_loc, lemur_data_loc, out, working, number_reads_generated):
    
    lemur_data = pd.read_csv(lemur_data_loc, delimiter='\t')
    print(lemur_data)
    metadata = pd.read_csv(metadata_loc)
    genome_list = metadata[metadata['Presence/Absence'] == 'Present']
    
    ## get the abundances from lemur and normalize to sum to 100
    genome_list['Abundance'] = genome_list['Taxonomy ID'].apply(lambda row: get_lemur_abundance(row, lemur_data))
    
    total_abundance = genome_list['Abundance'].sum()
    if total_abundance > 0:
        genome_list['Abundance'] = np.round(genome_list['Abundance'] / total_abundance * 100, 2)
    
    ## fix the locations names and drop unnecessary taxonomy column
    genome_list['Assembly Accession ID'] = f'{working}/magnet/reference_genomes/' + genome_list['Assembly Accession ID'].astype(str) + '.fasta'

    ## select final columns
    genome_list = genome_list[['Organism of Assembly', 'Assembly Accession ID', 'Abundance']]
    
    abundances =  genome_list[['Organism of Assembly', 'Abundance']] 
    
    out_loc = os.path.join(out, 'abundances.tsv')
    abundances.to_csv(out_loc, sep='\t', header=['Size', number_reads_generated], index=False)
    
    out_loc = os.path.join(out, 'genome_list1.tsv')
    out_loc2 = os.path.join(out, 'genome_list2.tsv')
    genome_list.to_csv(out_loc, sep='\t', header=False, index=False)
    genome_list[['Organism of Assembly', 'Assembly Accession ID']].to_csv(out_loc2, sep='\t', header=False, index=False)
    
    return genome_list
    
def make_genome_list(metadata_loc, out):
    metadata = pd.read_csv(metadata_loc)
    genome_list = metadata[['Organism of Assembly', 'Assembly Accession ID']]
    genome_list['Assembly Accession ID'] = 'test/magnet/reference_genomes/' + genome_list['Assembly Accession ID'].astype(str) + '.fasta'
    
    out_loc = os.path.join(out, 'genome_list.tsv')
    genome_list.to_csv(out_loc, sep='\t', header=False, index=False)
    
    return genome_list


def get_kraken_abundance(taxid, kraken_data):
    try:
        out = kraken_data[kraken_data['TaxID'] == taxid]
        if not out.empty:
            return out['Abundance'].values[0]
        else:
            return 0
    except KeyError:
        return 0

def get_final_species_abundances(kraken_report, metadata_loc, out, number_reads_generated):
    # Read kraken report and metadata
    kraken_data = pd.read_csv(kraken_report, delimiter='\t', header=None, names=['Abundance', 'NumCovered', 'NumTaxon', 'Rank', 'TaxID', 'Name'])
    metadata = pd.read_csv(metadata_loc)
    
    metadata = metadata[metadata['Presence/Absence'] == 'Present']
    metadata['Abundance'] = metadata['Taxonomy ID'].apply(lambda row: get_kraken_abundance(row, kraken_data))
        
    total_abundance = metadata['Abundance'].sum()
    if total_abundance > 0:
        metadata['Abundance'] = np.round(metadata['Abundance'] / total_abundance * 100, 2)

    abundances =  metadata[['Organism of Assembly', 'Abundance']] 
    
    out_loc = os.path.join(out, 'abundances.tsv')
    abundances.to_csv(out_loc, sep='\t', header=['Size', number_reads_generated], index=False)
      
    return abundances

def generate_species_file_info(genome_list, out):
    
    results = []
    for _, row in genome_list.iterrows():
        species = row['Organism of Assembly']
        filename = row['Assembly Accession ID']  
        
        if os.path.exists(filename):
            with open(filename, 'r') as file:
                first_line = file.readline().strip()[1:]
        else:
            first_line = 'File not found'
        
        results.append([species, first_line, 'circular'])
    
    species_info = pd.DataFrame(results, columns=['Species', 'FirstLine', 'Circular'])
    
    out_loc = os.path.join(out, 'species_info.tsv')
    species_info.to_csv(out_loc, sep='\t', header=False, index=False)
    
    return species_info


