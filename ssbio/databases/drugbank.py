from ssbio.utils import Date
date = Date()

import os
import os.path as op
import requests
import zipfile

drugbank_xml = 'http://www.drugbank.ca/releases/4-5-0/downloads/all-full-database'
drugbank_xsd = 'http://www.drugbank.ca/docs/drugbank.xsd'


def download_drugbank_db(username, password, outdir=os.getcwd()):
    # TODO: obfuscate password?
    r = requests.get(drugbank_xml, auth=(username, password))

    if r.headers['Content-Type'] == 'text/html; charset=utf-8':
        raise ConnectionError('DrugBank data could not be downloaded. Please check your username and password.')

    outfile = "{}-drugbank_all_full_database.xml.zip".format(date.short_date)
    if outdir:
        outfile = op.join(outdir, outfile)

    with open(outfile, "wb") as db:
        db.write(r.content)

    zip_ref = zipfile.ZipFile(outfile, 'r')
    zip_ref.extractall(outdir)
    zip_ref.close()

    return zip_ref


if __name__ == '__main__':
    import argparse

    p = argparse.ArgumentParser(description='DrugBank tools')
    p.add_argument('user', help='DrugBank user name')
    p.add_argument('passwd', help='DrugBank password')
    p.add_argument('--outdir', '-o', help='Directory to write to')
    args = p.parse_args()

    download_drugbank_db(args.user, args.passwd, args.outdir)