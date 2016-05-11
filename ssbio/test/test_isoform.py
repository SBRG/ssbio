import unittest

import ssbio.databases.isoform


class TestIsoform(unittest.TestCase):
    """Unit tests for PDB
    """

    def test_gene_isos_to_uniprot(self):
        gene_1 = '52'
        isos_1 = ['1','2','3']
        uniprot_isos_1 = ['P24666-1','P24666-2','P24666-3']
        unrev_uniprot_isos_1 = ['A0A0B4J207-1']
        expected_1 = {'1':'P24666-1','2':'P24666-2','3':'P24666-3'}
        self.assertEqual(ssbio.databases.isoform.gene_isos_to_uniprot(gene_1, isos_1, uniprot_isos_1, unrev_uniprot_isos_1), expected_1)

        gene_2 = '52'
        isos_2 = ['1', '2', '3']
        uniprot_isos_2 = ['P24666-1', 'P24666-2']
        unrev_uniprot_isos_2 = []
        expected_2 = {'1': 'P24666-1', '2': 'P24666-2', '3': None}
        self.assertEqual(ssbio.databases.isoform.gene_isos_to_uniprot(gene_2, isos_2, uniprot_isos_2, unrev_uniprot_isos_2), expected_2)

        gene_3 = '221823'
        isos_3 = ['1', '2']
        uniprot_isos_3 = ['P21108-1']
        unrev_uniprot_isos_3 = ['A0A0B4J207-1']
        expected_3 = {'1': 'P21108-1', '2': 'A0A0B4J207-1'}
        self.assertEqual(ssbio.databases.isoform.gene_isos_to_uniprot(gene_3, isos_3, uniprot_isos_3, unrev_uniprot_isos_3), expected_3)