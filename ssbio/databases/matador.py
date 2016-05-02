from ssbio import utils
date = utils.Date()

import os
import urllib.request
import io
import gzip
import os.path as op

matador_tsv_gz = "http://matador.embl.de/media/download/matador.tsv.gz"

def download_matador_db(outdir=None):
    if not outdir:
        outdir = os.getcwd()

    rename_unzipped_file = "{}-matador_full_database.tsv".format(date.short_date)

    if outdir:
        rename_unzipped_file = op.join(outdir, rename_unzipped_file)

    response = urllib.request.urlopen(matador_tsv_gz)
    compressed_file = io.BytesIO(response.read())
    decompressed_file = gzip.GzipFile(fileobj=compressed_file, mode='w')

    with open(rename_unzipped_file, "wt") as out_file:
        out_file.write(decompressed_file)

    return rename_unzipped_file