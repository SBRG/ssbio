import os.path as op
import tempfile
import unittest

import ssbio.databases.pdb

class TestPDB(unittest.TestCase):
    """Unit tests for PDB
    """

    def test_theoretical_pdbs(self):
        not_theoretical_pdbs = ['1kf6', '3bwm', '4ij5']
        theoretical_pdbs = ['1AAG','1AN3','1A8X','8CEL']

        # test if the non-theoretical pdbs return False when checked to be in the theoretical list
        for pdb in not_theoretical_pdbs:
            self.assertFalse(ssbio.databases.pdb.is_theoretical_pdb(pdb))

        # test if the theoretical pdbs return True when checked to be in the theoretical list
        for pdb in theoretical_pdbs:
            self.assertTrue(ssbio.databases.pdb.is_theoretical_pdb(pdb))

    def test_obsolete_pdbs(self):
        not_obsolete_pdbs = ['1kf6', '3bwm', '4ij5']
        obsolete_pdbs = ['5FDE', '4ZTW', '3V22', '1GEP']
        obsolete_pdbs_with_no_replacement = ['4EUR','4L2E','1FQH']

        # test if the non-obsolete pdbs return False when checked to be in the obsolete list
        for pdb in not_obsolete_pdbs:
            self.assertFalse(ssbio.databases.pdb.is_obsolete_pdb(pdb))
            self.assertIsNone(ssbio.databases.pdb.get_replacement_pdbs(pdb))

        # test if the obsolete pdbs return True when checked to be in the obsolete list
        for pdb in obsolete_pdbs:
            self.assertTrue(ssbio.databases.pdb.is_obsolete_pdb(pdb))
            self.assertIsInstance(ssbio.databases.pdb.get_replacement_pdbs(pdb), list)

        for pdb in obsolete_pdbs_with_no_replacement:
            self.assertTrue(ssbio.databases.pdb.is_obsolete_pdb(pdb))
            self.assertIsNone(ssbio.databases.pdb.get_replacement_pdbs(pdb))

    def test_sifts_pdb_chain_to_uniprot(self):
        mapping = {('104l', 'A'): ['P00720'],
                   ('10gs', 'A'): ['P09211']
                   }
        for k,v in mapping.items():
            self.assertEqual(ssbio.databases.pdb.sifts_pdb_chain_to_uniprot(k[0], k[1]), v)

    def test_top_pdb_blast_hit(self):
        seq = 'VLSPADKTNVKAAWGVKALSPADKTNVKAALTAVAHVDDMPNAL'
        evalue = 1
        results = {'hit_evalue': '0.0640608',
                     'hit_num_ident': '36',
                     'hit_pdb': '3IA3',
                     'hit_pdb_chain': 'B,D',
                     'hit_percent_ident': 0.81818181818181823,
                     'hit_score': '76'}
        self.assertEqual(ssbio.databases.pdb.top_pdb_blast_hit(seq,evalue=evalue), results)