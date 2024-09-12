import os, datetime

filemtime2datetime=lambda mypath: datetime.datetime.fromtimestamp(os.stat(mypath).st_mtime)

def ncbi_taxonomy_parse_nodes_dmp(nodes_dmp_path):
    '''
    Opens the NCBI nodes.dmp file from the taxdump FTP site. Parses it into a dictionary with
    key/vals in the form:
        { <taxon_id>: (<parent_taxon_id>, <level>), ...}

    Args:
        nodes_dmp_path (str): path to the file `nodes.dmp` from the NCBI taxonomy.

    Returns:
        :obj:`dict` object representing the nodes.dmp file. Key=NCBI Taxon ID, Val=tuple of (Taxon ID of Parent, Taxon Rank
                of self).
    '''
    if not os.path.isfile(nodes_dmp_path):
        print('The path specified in for the ncbi taxonomy file could not be opened.')
        return
    with open(nodes_dmp_path, 'r') as ncbi_f:
        # ncbi_dict = dict(map(lambda x: (int(x[0]), (int(x[2]), x[4], int(x[30]))), map(lambda x: x.strip().split('\t'), ncbi_f.readlines())))
        ncbi_dict = dict(map(lambda x: (int(x[0]), (int(x[2]), x[4])),
                             map(lambda x: x.strip().split('\t'), ncbi_f.readlines())))
    return ncbi_dict

def ncbi_taxonomy_parse_names_dmp(names_dmp_path):
    '''
    Opens the NCBI names.dmp file from the taxdump FTP site. Parses it into a dictionary with
    key/vals in the form: { <taxon_id>: (<taxon_name>, <name_alt>, <name_type>) }. Filters these
    results to only include lines that are of type "scientific name".

    Args:
        names_dmp_path (str): path to the file 'names.dmp' from the NCBI taxonomy

    Returns:
        dict: names_dict object as described above
    '''
    if not os.path.isfile(names_dmp_path):
        print(f'The path specified is not a valid file path: {names_dmp_path}')
        return {}
    with open(names_dmp_path,'r') as names_dmp_f:
        names_dmp = [i.strip().split('\t') for i in names_dmp_f.readlines()]
    return {int(x[0]): (x[2], x[4], x[6]) for x in names_dmp if x[6]=='scientific name'}

class NCBItaxonomy():
    '''Helper class to hold the NCBI taxonomy names/nodes files and easily return info without having to pass
    the big dictionary objects in every time.'''
    rank_main_seven = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    rank_prefixes = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

    def __init__(self, taxdmp_folder, default_dmp_filenames=True, lazy_parse=False):
        '''
        Initializes the object. Cannonically the taxdmp object from the ftp site is all in one folder, so this thing
        can be initialized just by specifying that.
        '''
        self.folder = taxdmp_folder
        if default_dmp_filenames:
            self.names_filepath = os.path.join(self.folder, 'names.dmp')
            self.nodes_filepath = os.path.join(self.folder, 'nodes.dmp')
            self.names_fetch_date = filemtime2datetime(self.names_filepath) if os.path.isfile(self.names_filepath) else None
            self.nodes_fetch_date = filemtime2datetime(self.nodes_filepath) if os.path.isfile(self.nodes_filepath) else None
            if not lazy_parse:
                self.parse_names_nodes_dmp()
            else:
                self._names_d = None
                self._nodes_d = None

    @property
    def nodes_d(self):
        if self._nodes_d is None:
            # logger.debug('self._nodes_d is None so importing now...')
            self._nodes_d = ncbi_taxonomy_parse_nodes_dmp(self.nodes_filepath)
        return self._nodes_d

    @property
    def names_d(self):
        if self._names_d is None:
            # logger.debug('self._names_d is None so importing now...')
            self._names_d = ncbi_taxonomy_parse_names_dmp(self.names_filepath)
        return self._names_d

    def parse_names_nodes_dmp(self):
        '''
        Just runs the external routines to parse the two files (provided they exist).
        '''
        # logger.info('Parsing both names/nodes files at object init.')
        self._names_d = ncbi_taxonomy_parse_names_dmp(self.names_filepath)
        self._nodes_d = ncbi_taxonomy_parse_nodes_dmp(self.nodes_filepath)

    def get_taxid_full_info(self, taxid, print_string=False):
        '''Returns the info from both files in a neatly organized format (i.e. TaxID, Name, Rank, Parent, Name-Type)'''
        nm_info = self.names_d[taxid]
        nd_info = self.nodes_d[taxid]
        full = {
            'taxid': taxid,
            'rank': nd_info[1],
            'name': nm_info[0],
            'parent': nd_info[0],
            'parent_name': self.names_d[nd_info[0]][0],
            'name-type': nm_info[-1]
        }
        pprint_output = '\n'.join([f'  {(k + ":"):<15s} {str(full[k]):s}' for k in full.keys()])
        if print_string:
            print(pprint_output)
        return full

    def get_taxid_lineage(self, taxid, return_only_main_ranks=False, format='list3tup'):
        '''Just a wrapper for the function below.'''
        L_three_tuples, main_7_rows = ncbi_taxonid_to_lineage_rawvector(taxid, self.nodes_d, self.names_d)
        if format=='list3tup':
            return L_three_tuples
        elif format=='rank2taxidnm':
            rank2taxidname = {i[1]: (i[0], i[2]) for i in L_three_tuples}
            return rank2taxidname

def ncbi_taxonid_to_lineage_rawvector(taxid, ncbi_taxonomy_nodes_dmp, ncbi_taxonomy_names_dmp=None):
    '''
    Takes a taxon-ID and returns a list of tuples containing the lineage from that taxon-ID up to the
    root. Returns a list sorted from highest-rank to lowest where each entry has the form:
        (Taxon-ID, Rank, Scientific-Name)

    Args:
        taxid (int):                Taxon-ID to look up.
        ncbi_taxonomy_nodes_dmp (dict): Dict object from parsing the ncbi taxonomy_nodes DEPRECATED
        return_only_main_ranks (bool):  If true, limits the vector to only the main 7 ranks (top=superkingdom now).

    Returns:
        lineage (list of tuples):   List of all the nodes in the tree in the form shown above.
    '''
    rank_main_seven = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    lineage = [taxid,]
    if taxid==0 or taxid not in ncbi_taxonomy_nodes_dmp:
        return []

    tdata = ncbi_taxonomy_nodes_dmp[taxid] # (parent_id, rank)
    ranks = [tdata[1],]      # ranks   = [taxid_rank, ]
    lineage.append(tdata[0]) # lineage = [taxid,     parent_taxid,]
    next = tdata[0]
    level_ct = 0
    while next > 1:
        tdata = ncbi_taxonomy_nodes_dmp[next]
        ranks.append(tdata[1])
        lineage.append(tdata[0])
        next = tdata[0]
        level_ct += 1
        if level_ct > 100:
            print ("taxid %s has a lineage that is supposedly 100+ levels")
            break
    if next==1:
        ranks.append(ncbi_taxonomy_nodes_dmp[next][1])
    #
    if ncbi_taxonomy_names_dmp is not None:
        taxon_names = ['' if ncbi_taxonomy_names_dmp is None else ncbi_taxonomy_names_dmp[i][0] for i in lineage]
    # Return it in descending order of the hierarchy:
    list3tup = list(zip(lineage, ranks, taxon_names))[::-1]
    main_rank_rows = [ranks[::-1].index(i) if i in ranks else -1 for i in rank_main_seven]
    return list3tup, main_rank_rows

