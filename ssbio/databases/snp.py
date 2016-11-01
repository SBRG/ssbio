import tempfile
import os.path as op
from requests import get  # to make GET request
import pandas as pd
import gzip
import urllib.request


def download(url, file_name):
    # open in binary mode
    with open(file_name, "wb") as file:
        # get request
        response = get(url)
        # write to file
        file.write(response.content)
    return file_name


def download_parse_uniprot_humsavar(wd):
    # TODO: this is not functional as of now - the whitespace separator doesn't work for the last column
    baseURL = "http://www.uniprot.org/docs/"
    filename = "humsavar.txt"

    downloaded_file = download(baseURL + filename, op.join(tempfile.gettempdir(), filename))
    parsed_file = pd.read_table(downloaded_file, skiprows=30, skipfooter=4, sep=r"\s*",
                                header=None, engine='python').rename(columns={0: 'u_gene_name', 1: 'u_uniprot_acc',
                                                                              2: 'u_ft_id', 3: 'u_var', 4: 'u_var_type',
                                                                              5: 'u_dbsnp', 6: 'u_disease'})
    outfile = op.join(wd, 'humsavar.csv')
    parsed_file.to_csv(outfile)
    return outfile


def download_parse_uniprot_variation(wd):
    baseURL = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/"
    filename = "homo_sapiens_variation.txt.gz"

    req = urllib.request.urlopen(baseURL + filename)
    decompressedFile = gzip.GzipFile(fileobj=req)

    outfile = op.join(tempfile.gettempdir(), filename)
    with open(outfile, 'w') as f:
        f.write(decompressedFile.read())

    uniprot_snps_new_df = pd.read_table(outfile, skiprows=139, delimiter='\t', header=None).rename(
        columns={0: 'u_gene_name', 1: 'u_uniprot_acc', 2: 'u_aa_change', 3: 'u_dbsnp',
                 4: 'u_var_type', 5: 'u_consq_type', 6: 'u_clin_sig', 7: 'u_phenotype',
                 8: 'u_cytogenic_band', 9: 'u_chr_loc', 10: 'ensg', 11: 'enst', 12: 'ensp'})

    outfile_pandas = op.join(wd, 'homo_sapiens_variation.csv')
    uniprot_snps_new_df.to_csv(outfile_pandas)
    return outfile_pandas


if __name__ == "__main__":
    pass