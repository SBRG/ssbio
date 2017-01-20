import unittest

import ssbio.sequence.utils.alignment as ssbaln
from Bio.Align import MultipleSeqAlignment

class TestAlignment(unittest.TestCase):
    """Unit tests for Alignment"""

    def test_pairwise_sequence_alignment_str(self):
        base_id = 'a_test'
        base = 'MSPRVGVTLSGRYRLQRLIATGGMGQVWEAVDNRLGRRVAVVASASA'
        muts_id = 'b_test'
        muts = 'MSPRVGVTLSGRYRLQRLIATGGMGQVWEAVDNRLGRRVARRVAASSSAS'

        alignment = ssbaln.pairwise_sequence_alignment(a_seq_id=base_id, a_seq=base,
                                                       b_seq_id=muts_id, b_seq=muts,
                                                       engine='biopython')
        self.assertTrue(isinstance(alignment, MultipleSeqAlignment))
        self.assertEqual(alignment[0].id, base_id)
        self.assertEqual(alignment[0].seq, 'MSPRVGVTLSGRYRLQRLIATGGMGQVWEAVDNRLGRRVAV--VA-SASA-')
        self.assertEqual(alignment[1].id, muts_id)
        self.assertEqual(alignment[1].seq, 'MSPRVGVTLSGRYRLQRLIATGGMGQVWEAVDNRLGRRVA-RRVAASSSAS')
        self.assertEqual(alignment.annotations, {'end': 51, 'percent_identity': 88.23529411764706, 'score': 224.0, 'start': 0})

        alignment2 = ssbaln.pairwise_sequence_alignment(a_seq_id=base_id, a_seq=base,
                                                        b_seq_id=muts_id, b_seq=muts,
                                                        engine='needle', force_rerun=True)
        self.assertTrue(isinstance(alignment2, MultipleSeqAlignment))
        self.assertEqual(alignment2[0].id, base_id)
        self.assertEqual(alignment2[0].seq, 'MSPRVGVTLSGRYRLQRLIATGGMGQVWEAVDNRLGRRVA--VVASASA-')
        self.assertEqual(alignment2[1].id, muts_id)
        self.assertEqual(alignment2[1].seq, 'MSPRVGVTLSGRYRLQRLIATGGMGQVWEAVDNRLGRRVARRVAASSSAS')
        self.assertEqual(alignment2.annotations, {'score': 213.5, 'percent_gaps': 6.0, 'percent_similarity': 92.0, 'percent_identity': 90.0})