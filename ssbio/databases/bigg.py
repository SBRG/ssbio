import json
import requests
import tempfile
import ssbio.databases.pdb
import ssbio.utils
import os.path as op
from bioservices.uniprot import UniProt
import ssbio.databases.uniprot
bs_unip = UniProt()


def get_pdbs_for_gene(bigg_model, bigg_gene, cache_dir=tempfile.gettempdir(), force_rerun=False):
    """Attempt to get a rank-ordered list of available PDB structures for a BiGG Model and its gene.

    Args:
        bigg_model: BiGG Model ID
        bigg_gene: BiGG Gene ID

    Returns:
        list: rank-ordered list of tuples of (pdb_id, chain_id)

    """
    my_structures = []

    # Download gene info
    gene = ssbio.utils.request_json(link='http://bigg.ucsd.edu/api/v2/models/{}/genes/{}'.format(bigg_model, bigg_gene),
                                    outfile='{}_{}.json'.format(bigg_model, bigg_gene),
                                    outdir=cache_dir,
                                    force_rerun_flag=force_rerun)

    uniprots = []
    if 'database_links' in gene:
        if 'UniProt' in gene['database_links']:
            uniprots = [x['id'] for x in gene['database_links']['UniProt']]
        elif 'NCBI GI' in gene['database_links']:
            uniprots = []
            gis = [x['id'] for x in gene['database_links']['NCBI GI']]
            gi_uniprots = bs_unip.mapping(fr='P_GI', to='ACC', query=gis).values()
            uniprots.extend(gi_uniprots)
            uniprots = ssbio.utils.flatlist_dropdup(uniprots)
            uniprots = [x for x in uniprots if ssbio.databases.uniprot.is_valid_uniprot_id(x)]

    if uniprots:
        for u in uniprots:
            get_best_structure = ssbio.databases.pdb.best_structures(uniprot_id=u, outdir=cache_dir)
            if get_best_structure:
                for best_structure in get_best_structure:
                    my_structures.append((best_structure['pdb_id'], best_structure['chain_id']))

    return my_structures
