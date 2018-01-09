"""
SWISSMODEL
==========
"""

import json
import logging
import requests
import ssbio.utils
import os.path as op
from collections import defaultdict

log = logging.getLogger(__name__)


class SWISSMODEL():
    """Methods to parse through a SWISS-MODEL metadata set.

    Download a particular organism's metadata from SWISS-MODEL here: https://swissmodel.expasy.org/repository

    Args:
        metadata_dir (str): Path to the extracted SWISS-MODEL_Repository folder

    """

    def __init__(self, metadata_dir):
        self.metadata_dir = metadata_dir
        """str: Path to the extracted SWISS-MODEL_Repository folder"""

        self.all_models = None
        """dict: Dictionary of lists, UniProt ID as the keys"""

        # Parse the INDEX_JSON file and then store all the metadata in all_models
        self.parse_metadata()

    @property
    def metadata_index_json(self):
        """str: Path to the INDEX_JSON file."""
        return op.join(self.metadata_dir, 'INDEX_JSON')

    @property
    def uniprots_modeled(self):
        """list: Return all UniProt accession numbers with at least one model"""
        return list(self.all_models.keys())

    def parse_metadata(self):
        """Parse the INDEX_JSON file and reorganize it as a dictionary of lists."""

        all_models = defaultdict(list)

        with open(self.metadata_index_json) as f:
            loaded = json.load(f)

        for m in loaded:
            all_models[m['uniprot_ac']].append(m)

        self.all_models = dict(all_models)

    def get_models(self, uniprot_acc):
        """Return all available models for a UniProt accession number.

        Args:
            uniprot_acc (str): UniProt ACC/ID

        Returns:
            dict: All available models in SWISS-MODEL for this UniProt entry

        """
        return self.all_models[uniprot_acc]

    def download_models(self, uniprot_acc, outdir='', force_rerun=False):
        """Download all models available for a UniProt accession number.

        Args:
            uniprot_acc (str): UniProt ACC/ID
            outdir (str): Path to output directory, uses working directory if not set
            force_rerun (bool): Force a redownload the models if they already exist

        Returns:
            list: Paths to the downloaded models

        """
        downloaded = []
        subset = self.get_models(uniprot_acc)

        for entry in subset:
            ident = '{}_{}_{}_{}'.format(uniprot_acc, entry['template'], entry['from'], entry['to'])
            outfile = op.join(outdir, ident + '.pdb')

            if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
                response = requests.get(entry['url'])

                if response.status_code == 404:
                    log.error('{}: 404 returned, no model available.'.format(ident))

                else:
                    with open(outfile, 'w') as f:
                        f.write(response.text)

                    log.debug('{}: downloaded homology model'.format(ident))
                    downloaded.append(outfile)
            else:
                downloaded.append(outfile)

        return downloaded