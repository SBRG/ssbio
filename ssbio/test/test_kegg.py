import unittest
import six

import ssbio.databases.kegg
from ssbio.databases.kegg import KEGGProp

class TestKEGG(unittest.TestCase):
    """Unit tests for KEGG
    """

    def test_kegg_prop(self):
        kg = KEGGProp(kegg_id='mtu:Rv0417',
                      sequence_file='test_files/sequences/mtu-Rv0417.faa',
                      metadata_file='test_files/sequences/mtu-Rv0417.kegg')

        kg_keys = ['bigg','gene_name','id','kegg','metadata_file','metadata_path','pdbs','refseq',
                   'sequence_file','seq_len','sequence_path','taxonomy','uniprot']

        for k in kg_keys:
            self.assertTrue(hasattr(kg, k))

        self.assertEqual(kg.bigg, None)
        self.assertEqual(kg.gene_name, None)
        self.assertEqual(kg.id, 'mtu:Rv0417')
        self.assertEqual(kg.kegg, 'mtu:Rv0417')
        self.assertEqual(kg.metadata_file, 'mtu-Rv0417.kegg')
        self.assertEqual(kg.pdbs, None)
        self.assertEqual(kg.refseq, 'NP_214931')
        self.assertEqual(kg.sequence_file, 'mtu-Rv0417.faa')
        self.assertEqual(kg.seq_len, 252)
        self.assertEqual(kg.taxonomy, 'mtu  Mycobacterium tuberculosis H37Rv')
        six.assertCountEqual(self, kg.uniprot, ['P9WG73', 'I6Y3Q4'])