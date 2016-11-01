import requests


def chembl_uniprot_to_chembl_id(uniprot_id):
    """Map a UniProt ID to a ChEMBL ID.

    Args:
        uniprot_id (str): UniProt ID

    Returns:
        str: ChEMBL ID

    """
    r = requests.get('https://www.ebi.ac.uk/chemblws/targets/uniprot/{}.json'.format(uniprot_id))
    if r.status_code == 200:
        d = r.json()
        return str(d['target']['chemblId'])
    else:
        return None


def chembl_targets_to_approved_drugs(chembl_id):
    """Map a ChEMBL ID to the ChEMBL drugs that target it

    Args:
        chembl_id (str):

    Returns:
        ?

    """
    r = requests.get('https://www.ebi.ac.uk/chemblws/targets/{}/approvedDrug.json'.format(chembl_id))
    if r.status_code == 200:
        d = r.json()
        approved_drugs = d['approvedDrugs']
        if approved_drugs:
            drugs_names_actions = {}
            for drug in approved_drugs:
                drug_id = str(drug['chemblId'])
                drug_name = str(drug['name'])
                drug_mech = str(drug['mechanismOfAction'])
                drugs_names_actions[(drug_id,drug_name)] = drug_mech
            return drugs_names_actions
        else:
            return None


def chembl_uniprot_to_drugs(uniprot_id):
    """Map a UniProt ID to approved drugs that target it

    Args:
        uniprot_id (str): UniProt ID

    Returns:
        ?

    """
    chembl_id = chembl_uniprot_to_chembl_id(uniprot_id)
    if chembl_id:
        return chembl_targets_to_approved_drugs(chembl_id)
    else:
        return None