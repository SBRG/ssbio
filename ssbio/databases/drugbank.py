from ssbio import utils
date = utils.Date()

import os
import os.path as op
import pandas as pd
import xmltodict
import requests
import zipfile
import shutil
import collections
from tqdm import tqdm

drugbank_xml = 'http://www.drugbank.ca/releases/4-5-0/downloads/all-full-database'
drugbank_xsd = 'http://www.drugbank.ca/docs/drugbank.xsd'


# TODO: check out https://github.com/dhimmel/drugbank

def download_drugbank_db(username, password, outdir=None):
    """Download the DrugBank XML database and return the unzipped file path

    You need to have a DrugBank account now to download the XML file.
    Args:
        username: your DrugBank.ca email username
        password: your DrugBank.ca password
        outdir: optional output directory (default is the current directory)

    Returns:
        rename_unzipped_file (str): path to unzipped XML file

    """
    # TODO: obfuscate password?

    if not outdir:
        outdir = os.getcwd()

    # define filenames here
    outfile = "{}-drugbank_all_full_database.xml.zip".format(date.short_date)
    unzipped_file = 'full database.xml'
    rename_unzipped_file = "{}-drugbank_all_full_database.xml".format(date.short_date)

    if outdir:
        outfile = op.join(outdir, outfile)
        unzipped_file = op.join(outdir, unzipped_file)
        rename_unzipped_file = op.join(outdir, rename_unzipped_file)

    # if the file was already downloaded today, just return it
    if op.exists(rename_unzipped_file):
        return rename_unzipped_file

    r = requests.get(drugbank_xml, auth=(username, password))

    # an HTML page is returned if there is an error getting the file
    if r.headers['Content-Type'] == 'text/html; charset=utf-8':
        raise ConnectionError('DrugBank data_dir could not be downloaded. Please check your username and password.')

    with open(outfile, "wb") as db:
        db.write(r.content)

    zip_ref = zipfile.ZipFile(outfile, 'r')
    zip_ref.extractall(outdir)
    zip_ref.close()

    shutil.move(unzipped_file, rename_unzipped_file)

    return rename_unzipped_file


def parse_drugbank_xml(db_xml_path, outdir=None):
    """Parse the main DrugBank XML database

    Args:
        db_xml_path:
        outdir:

    Returns:
        Paths to Pandas DataFrames for:
        db_main_path: Main DrugBank ID mapping
        db_bind_path: DrugBank binding partners
        db_snp_path: Pharmacogenomics information

    TODO:
        - Use ElementTree instead

    """
    if not outdir:
        outdir = os.getcwd()

    db_main_path = op.join(outdir, '{}-db_main_df.csv'.format(date.short_date))
    db_bind_path = op.join(outdir, '{}-db_binding_df.csv'.format(date.short_date))
    db_snp_path = op.join(outdir, '{}-db_snp_df.csv'.format(date.short_date))

    if op.exists(db_main_path) and op.exists(db_bind_path) and op.exists(db_snp_path):
        return db_main_path, db_bind_path, db_snp_path

    db_xml = op.join(db_xml_path)
    db_xml_open = open(db_xml, "rb")
    db_dict = xmltodict.parse(db_xml_open)

    db_main = []
    db_bind_to = []
    db_snp = []

    for drug in tqdm(db_dict['drugbank']['drug']):  # [:100]:
        main_row = collections.OrderedDict()

        db_id = utils.force_list(drug['drugbank-id'])[0]['#text']

        main_row['db_id'] = db_id
        main_row['db_name'] = drug['name']
        main_row['db_type'] = drug['@type']
        main_row['db_groups'] = utils.force_string(drug['groups']['group'])

        # Percentage of the drug that is bound in plasma proteins
        main_row['db_protein_binding'] = drug['protein-binding']

        # PharmGKB ID
        if drug['external-identifiers']:
            for ext_id in utils.force_list(drug['external-identifiers']['external-identifier']):
                if ext_id['resource'] == 'PharmGKB':
                    main_row['db_pgkb'] = ext_id['identifier']

                    # these "reactions" are what the drug is metabolized into
                    #     print(row['db_id, drug['reactions'])

        db_main.append(main_row)

        # SNP-FX contains data on the drug, the interacting protein(kegg),
        # the ‘causal’ SNPs or genetic variants for that gene/protein,
        # the therapeutic response or effects caused by the SNP-drug interaction
        # (improved or diminished response, changed dosing requirements, etc.)
        # and the associated references describing these effects in more detail.
        # SNP-ADR follows a similar format to SNP-FX but the clinical responses
        # are restricted only to adverse drug reactions (ADR).
        # SNP-FX contains literature-derived data on the therapeutic effects or therapeutic responses
        # for more than 60 drug-polymorphism combinations,
        # while SNP-ADR contains data on adverse reactions compiled from more than 50 drug-polymorphsim pairings.
        # All of the data_dir in these tables is hyperlinked to drug entries from DrugBank, protein data_dir from UniProt,
        # SNP data_dir from dbSNP and bibliographic data_dir from PubMed.
        # NOTE: these dataframes are combined, with a db_snp_adr column indicating a known ADR
        if drug['snp-effects']:
            for snp_eff in utils.force_list(drug['snp-effects']['effect']):
                snp_row = collections.OrderedDict()

                snp_row['db_id'] = db_id
                snp_row['db_snp_uniprot'] = snp_eff['uniprot-id']
                snp_row['db_snp_rsid'] = snp_eff['rs-id']
                snp_row['db_snp_allele'] = snp_eff['allele']
                snp_row['db_snp_allele_info'] = snp_eff['defining-change']
                snp_row['db_snp_effect'] = snp_eff['description']
                snp_row['db_snp_adr'] = False

                db_snp.append(snp_row)

        if drug['snp-adverse-drug-reactions']:
            for snp_eff in utils.force_list(drug['snp-adverse-drug-reactions']['reaction']):
                snp_row = collections.OrderedDict()

                snp_row['db_id'] = db_id
                snp_row['db_snp_uniprot'] = snp_eff['uniprot-id']
                snp_row['db_snp_rsid'] = snp_eff['rs-id']
                snp_row['db_snp_allele'] = snp_eff['allele']
                snp_row['db_snp_allele_info'] = snp_eff['adverse-reaction']
                snp_row['db_snp_effect'] = snp_eff['description']
                snp_row['db_snp_adr'] = True

                db_snp.append(snp_row)

        for bind_type in ['targets', 'enzymes', 'carriers', 'transporters']:
            if drug[bind_type]:
                for bind_to in utils.force_list(drug[bind_type][bind_type.strip('kegg')]):
                    if 'polypeptide' in bind_to.keys():
                        for polypeptide in utils.force_list(bind_to['polypeptide']):
                            bind_to_row = collections.OrderedDict()
                            bind_to_row['db_id'] = db_id
                            bind_to_row['db_uniprot'] = polypeptide['@id']
                            bind_to_row['db_type'] = bind_type.strip('kegg')
                            bind_to_row['db_known_action'] = bind_to['known-action']
                            bind_to_row['db_location'] = polypeptide['cellular-location']
                            if bind_to['actions']:
                                bind_to_row['db_action'] = utils.force_string(bind_to['actions']['action'])

                            db_bind_to.append(bind_to_row)

    pd.DataFrame(db_main).to_csv(db_main_path)
    pd.DataFrame(db_bind_to).to_csv(db_bind_path)
    pd.DataFrame(db_snp).to_csv(db_snp_path)

    print('Saved parsed dataframes: {}, {}, {}'.format(db_main_path, db_bind_path, db_snp_path))
    return db_main_path, db_bind_path, db_snp_path

if __name__ == '__main__':
    import argparse

    p = argparse.ArgumentParser(description='DrugBank tools')
    p.add_argument('user', help='DrugBank user name')
    p.add_argument('passwd', help='DrugBank password')
    p.add_argument('--outdir', '-o', help='Directory to write to')
    args = p.parse_args()

    print('Downloading DrugBank XML...')
    db_xml_path = download_drugbank_db(args.user, args.passwd, args.outdir)
    print('Parsing DrugBank database into dataframes...')
    parse_drugbank_xml(db_xml_path, outdir=args.outdir)