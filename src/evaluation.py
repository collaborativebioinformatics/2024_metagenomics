
from src.tax_identification import *
from src.ncbi_taxonomy_utils.py import NCBItaxonomy

k2db='/home/Users/rdd4/VIMERA_DB_1.0/kraken2-genbank'
#working='/home/Users/mn56/work/mimic/results_1000_perf'
working='/home/Users/mn56/work/mimic/results_50000_perf'
#infasta='/home/Users/rdd4/Mimic/SRR30413550/simulated_data_SRR30413550/simulated1000.fasta'
#infasta='/home/Users/rdd4/Mimic/SRR30413550/simulated_data_SRR30413550/simulated1000_perf.fasta'
infasta='/home/Users/rdd4/Mimic/SRR30413550/simulated_data_SRR30413550/simulated50000_perf.fasta'
#ncbi_taxdmp_fold='/home/mnute/work/rice/ncbi_taxonomy'
ncbi_taxdmp_fold='/home/Users/mn56/data/ncbi_taxonomy'
magfold='/home/Users/rdd4/Mimic/SRR30413550/magnet'
k2_output_path='/home/Users/mn56/work/mimic/results_1000_perf/kraken2/output.txt'
k2_output_path_aln='/home/Users/mn56/work/mimic/results_1000/kraken2/output.txt'
k2_output_path_50k='/home/Users/mn56/work/mimic/results_50000_perf/kraken2/output.txt'
k2_output_path_50k_aln='/home/Users/mn56/work/mimic/results_50000/kraken2/output.txt'


run_kraken2(infasta, k2db, working)

k2score_50k=score_kraken2_nanosim_output(k2_output_path_50k, magfold, ncbi_taxdmp_fold)
rescore_kraken2_nanosim_output_by_rank(k2score_50k[3], ncbi_taxdmp_fold)

def get_prec_rec(fn_tp_fp_tuple):
    return fn_tp_fp_tuple[1]/(fn_tp_fp_tuple[1]+fn_tp_fp_tuple[2]), fn_tp_fp_tuple[1]/sum(fn_tp_fp_tuple[:3])

def read_magnet_cluster_representatives(magnet_folder):
    '''reads the csv file 'cluster_representative.csv' within the magnet subfolder'''
    filepath=os.path.join(magnet_folder,'cluster_representative.csv')
    lines=[]
    with open(filepath,'r') as magcr:
        lines=magcr.readlines()
        headers=lines[0].strip().split()
        clustreps=[i.strip().split(',') for i in lines[1:]]
        clust_id_lkp = {i[5]: int(i[0]) for i in clustreps}
        #species_id_lkp = {i[9]: int(i[0]) for i in clustreps}
    return clust_id_lkp


def score_kraken2_nanosim_output(kraken2_output, magnet_folder):
    '''parses a standard kraken2 output file and reports results. assuming for now that
    all classifications and all ground truth taxon IDs are at the species level.'''
    # Get kraken2 output.txt file
    with open(kraken2_output,'r') as k2out:
        k2output = [i.strip().split('\t') for i in k2out.readlines()]
    
    # Get magnet cluster_represnatives taxon ID lookup:
    magcr_lkp = read_magnet_cluster_representatives(magnet_folder)
    #
    res = []
    for i in k2output:
        read_hdr = i[1]
        is_classified = i[0]
        taxon_id = int(i[2])
        ground_truth_organism = read_hdr.split('-')[0].replace('_',' ')
        true_taxon_id = magcr_lkp[ground_truth_organism]
        res.append((read_hdr, is_classified, true_taxon_id, taxon_id))
    
    false_negatives=sum([1 if i[1]=='U' else 0 for i in res])
    true_positives=sum([1 if i[2]==i[3] else 0 for i in res])
    false_positives = sum([1 if i[1]=='C' and i[2]!=i[3] else 0 for i in res])
    return false_negatives, true_positives, false_positives, res

ncbitax=NCBItaxonomy(ncbi_taxdmp_fold)

def rescore_kraken2_nanosim_output_by_rank(kraken2_scores, ncbitax):
    '''takes the list of tuples from the previous function where indices 2,3 
    of the tuple are the true taxon ID and the kraken estimated taxon ID, 
    respectively, and restates the FN/TP/FP counts by rank.'''
    #
    def fn_tp_fp_tn(id_pair_list):
        tn = sum([1 if i[0] is None and i[1] is not None else 0 for i in id_pair_list])
        fn = sum([1 if i[1] is None and i[0] is not None else 0 for i in id_pair_list])
        tp = sum([1 if i[0] is not None and i[1] is not None and i[0]==i[1] else 0 for i in id_pair_list])
        fp = sum([1 if i[0] is not None and i[1] is not None and i[0]!=i[1] else 0 for i in id_pair_list])
        prec = tp / (tp+fp)
        rec = tp / (tp+fp+fn)
        return fn, tp, fp, tn, prec, rec
    #
    species_res=[]
    genus_res = []
    family_res = []
    order_res = []
    class_res = []
    phylum_res = []
    for i in kraken2_scores:
        true_taxid = i[2]; est_taxid=i[3];
        try:
            true_lineage = ncbitax.get_taxid_lineage(true_taxid, True, 'rank2taxidnm')
        except:
            true_lineage = {}
        try:
            if est_taxid==0:
                est_lineage={}
            else:
                est_lineage = ncbitax.get_taxid_lineage(est_taxid, True, 'rank2taxidnm')
        except:
            est_lineage={}
        species_res.append((true_lineage.get('species',None), est_lineage.get('species',None)))
        genus_res.append((true_lineage.get('genus',None), est_lineage.get('genus',None)))
        family_res.append((true_lineage.get('family',None), est_lineage.get('family',None)))
        order_res.append((true_lineage.get('order',None), est_lineage.get('order',None)))
        class_res.append((true_lineage.get('class',None), est_lineage.get('class',None)))
        phylum_res.append((true_lineage.get('phylum',None), est_lineage.get('phylum',None)))
    #
    full_res = {}
    full_res['species'] = fn_tp_fp_tn(species_res)
    full_res['genus'] = fn_tp_fp_tn(genus_res)
    full_res['family'] = fn_tp_fp_tn(family_res)
    full_res['order'] = fn_tp_fp_tn(order_res)
    full_res['class'] = fn_tp_fp_tn(class_res)
    full_res['phylum'] = fn_tp_fp_tn(phylum_res)
    #
    print('rank\tFN\tTP\tFP\tTN\tPrec\tRec')
    for i in full_res.keys():
        print('%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f' % ((i,)+full_res[i]))
    #
    return full_res

