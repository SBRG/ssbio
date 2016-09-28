import pandas as pd
import requests
from bs4 import BeautifulSoup
import numpy as np

# this takes a uniprot id returns its chembl id
def chembl_uniprot_to_chembl_id(uniprot_id):
    r = requests.get('https://www.ebi.ac.uk/chemblws/targets/uniprot/%s.json' % uniprot_id)
    if r.status_code == 200:
        d = r.json()
        return str(d['target']['chemblId'])
    else:
        return None


# this takes a chembl ID and returns drugs that target it
def chembl_targets_to_drugs(chembl_id):
    r = requests.get('https://www.ebi.ac.uk/chemblws/targets/%s/approvedDrug.json' % chembl_id)
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


# wrapper function for the above 2
def chembl_uniprot_to_drugs(uniprot_id):
    chembl_id = chembl_uniprot_to_chembl_id(uniprot_id)
    if chembl_id:
        return chembl_targets_to_drugs(chembl_id)
    else:
        return None