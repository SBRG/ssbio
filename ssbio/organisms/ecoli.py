import math
import warnings

import cachetools
from tqdm import tqdm

import ssbio.databases.identifiers
import ssbio.databases.uniprot
import ssbio.sequence.properties.aggregation_propensity as agg
import ssbio.sequence.properties.kinetic_folding_rate as kfr
import ssbio.sequence.properties.thermostability as ts


def convert_bnumber_to_uniprot(bnumber):
    """Map an E. coli locus id (ie. b0003) to its reviewed UniProt ID

    Args:
        bnumber: E. coli locus ID

    Returns:
        UniProt ID (reviewed only) - in most cases there should only be one
    """
    genes_to_uniprots = ssbio.databases.identifiers.bioservices_uniprot_mapper('ENSEMBLGENOME_ID', 'ACC', bnumber)
    reviewed_uniprots = [x for x in genes_to_uniprots[bnumber] if ssbio.databases.uniprot.uniprot_reviewed_checker(x)]
    if len(reviewed_uniprots) > 1:
        warnings.warn('{}: more than one reviewed UniProt entry. Returning first mapped ID.'.format(bnumber))
    return reviewed_uniprots[0]

def convert_bnumbers_to_uniprot_batch(bnumbers):
    """Map a list of E. coli locus ids (ie. b0003) to their reviewed UniProt IDs

    Args:
        bnumbers: list of E. coli locus IDs

    Returns:
        Dictionary of locus ID to UniProt ID (reviewed only) - in most cases there should only be one
    """
    final_mapping = {}
    genes_to_uniprots = ssbio.databases.identifiers.bioservices_uniprot_mapper('ENSEMBLGENOME_ID', 'ACC', bnumbers)
    for k,v in genes_to_uniprots.items():
        reviewed_uniprots = [x for x in v if ssbio.databases.uniprot.uniprot_reviewed_checker(x)]
        if len(reviewed_uniprots) > 1:
            warnings.warn('{}: more than one reviewed UniProt entry. Returning first mapped ID.'.format(k))
        final_mapping[k] = reviewed_uniprots[0]
    return final_mapping

@cachetools.func.ttl_cache(maxsize=500)
def uniprot_info(bnumber):
    """Get the mapped UniProt ID and sequence for an E. coli locus id

    Args:
        bnumber: E. coli locus ID

    Returns:
        tuple of (uniprot ID, sequence)
    """
    uniprot_id = convert_bnumber_to_uniprot(bnumber)
    uniprot_seq = ssbio.databases.uniprot.get_fasta(uniprot_id)
    return(uniprot_id, uniprot_seq)

def get_dG_for_gene(bnumber, temp):
    """Predict dG at temperature T for an E. coli b number

    1. Maps the b-number to the corresponding UniProt ID
    2. Retrieves the sequence for this ID
    3. Calculates the dG and Keq

    Args:
        bnumber: E. coli locus ID
        temp: temperature in celcius

    Returns:
        (dG (cal/mol), Keq, method)

    """
    id_and_seq = uniprot_info(bnumber)
    return ts.get_dG_at_T(id_and_seq[1], temp)

def get_dG_for_genes_batch(bnumbers, temp):
    """Predict dG at temperature T for a list of E. coli b numbers

    1. Maps the b-number to the corresponding UniProt ID
    2. Retrieves the sequence for this ID
    3. Calculates the dG and Keq

    Args:
        bnumbers: list of E. coli locus IDs
        temp: temperature in celcius

    Returns:
        {bnumber: (dG (cal/mol), Keq, method)}

    """
    final_dict = {}
    id_mapping = convert_bnumbers_to_uniprot_batch(bnumbers)
    for bnumber,uniprot in tqdm(id_mapping.items()):
        id_and_seq = uniprot_info(bnumber)
        final_dict[bnumber] = ts.get_dG_at_T(id_and_seq[1], temp)
    return final_dict


@cachetools.func.ttl_cache(maxsize=500)
def get_aggregation_propensity_for_gene(bnumber):
    """Predict aggregation propensity

    Args:
        bnumber: E. coli locus ID

    Returns:
        ?

    """
    agg_predictions = agg.AMYLPRED('kec003@ucsd.edu', 'KoKa123456')
    id_and_seq = uniprot_info(bnumber)
    return agg_predictions.get_aggregation_propensity(id_and_seq[1])


@cachetools.func.ttl_cache(maxsize=500)
def get_folding_rate_for_gene(bnumber, temp, refT=37):
    """Predict the kinetic folding rate of a protein at temperature T

    Args:
        bnumber:  E. coli locus ID
        temp: temperature
        refT: reference temperature default to 37C

    Returns:
        kinetic folding rate kf (s-1)

    """

    slope = 22000  # not many data on this value, but it's effect on growth rate very small
    id_and_seq = uniprot_info(bnumber)

    # TODO: get (pdb_file, dssp_file, chain_id) for bnumber
    # GEMPRO
    # opn = ssbdssp.get_ss_class(pdb_file, dssp_file, chain)
    opn='mixed'

    # get folding rate for the reference temperature
    ref_rate = kfr.get_folding_rate_from_url(id_and_seq[1], opn)
    preFactor = float(ref_rate) + slope / (float(refT) + 273.15)

    # calculate folding rate at desired temperature
    rate = math.exp(preFactor - slope / (float(temp) + 273.15))

    return rate
