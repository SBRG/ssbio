import os.path as op
import tempfile
import unittest

import ssbio.databases.pdb

class TestPDB(unittest.TestCase):
    """Unit tests for PDB
    """
    # # TODO: this result changes based on pdb database
    # def test_top_pdb_blast_hit(self):
    #     seq = 'VLSPADKTNVKAAWGVKALSPADKTNVKAALTAVAHVDDMPNAL'
    #     evalue = 1
    #     results = {'hit_evalue': '0.0640608',
    #                  'hit_num_ident': '36',
    #                  'hit_pdb': '3IA3',
    #                  'hit_pdb_chain': 'B,D',
    #                  'hit_percent_ident': 0.81818181818181823,
    #                  'hit_score': '76'}
    #     blast_hit = ssbio.databases.pdb.top_pdb_blast_hit(seq, evalue=evalue)
    #     for x in results.keys():
    #         self.assertTrue(x in list(blast_hit.keys()))