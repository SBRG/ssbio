import unittest

import ssbio.sequence.properties.thermostability


class TestThermostability(unittest.TestCase):
    """Unit tests for thermostability calculations
    """

    def test_gene_isos_to_uniprot(self):
        gene_1 = '52'
        isos_1 = ['1', '2', '3']
        uniprot_isos_1 = ['P24666-1', 'P24666-2', 'P24666-3']
        unrev_uniprot_isos_1 = ['A0A0B4J207-1']
        expected_1 = {'1': 'P24666-1', '2': 'P24666-2', '3': 'P24666-3'}
        self.assertEqual(
            ssbio.databases.isoform.gene_isos_to_uniprot(gene_1, isos_1, uniprot_isos_1, unrev_uniprot_isos_1),
            expected_1)