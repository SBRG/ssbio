import pytest
from ssbio.core.protein import Protein
import os.path as op


@pytest.fixture(scope='class')
def test2(test_files_sequences, test_files_structures, test_files_outputs):
    """Protein with confusing mappings between residue numbers"""
    p = Protein(ident='P21964', root_dir=test_files_outputs)
    p.load_uniprot(uniprot_id='P21964',
                   uniprot_seq_file=op.join(test_files_sequences, 'P21964.fasta'),
                   uniprot_xml_file=op.join(test_files_sequences, 'P21964.xml'),
                   download=False, set_as_representative=True)
    p.load_pdb(pdb_id='3bwm', mapped_chains=['A'], pdb_file=op.join(test_files_structures, '3bwm.pdb'),
               file_type='pdb', representative_chain='A', set_as_representative=True)
    p.align_seqprop_to_structprop(seqprop=p.sequences.get_by_id('P21964'),
                                  structprop=p.structures.get_by_id('3bwm'),
                                  chain_id='A')
    return p



@pytest.fixture(scope='class')
def straightforward_resnum_mapping(test_files_sequences, test_files_structures, test_files_outputs):
    """Protein with confusing mappings between residue numbers"""
    p = Protein(ident='P0ABP8', root_dir=test_files_outputs)
    p.load_uniprot(uniprot_id='P0ABP8',
                   uniprot_seq_file=op.join(test_files_sequences, 'P0ABP8.fasta'),
                   uniprot_xml_file=op.join(test_files_sequences, 'P0ABP8.xml'),
                   download=False, set_as_representative=True)
    p.load_pdb(pdb_id='1a69', mapped_chains=['A'], pdb_file=op.join(test_files_structures, '1a69.pdb'),
               file_type='pdb', representative_chain='A', set_as_representative=True)
    p.align_seqprop_to_structprop(seqprop=p.sequences.get_by_id('P0ABP8'),
                                  structprop=p.structures.get_by_id('1a69'),
                                  chain_id='A')
    return p


@pytest.fixture(scope='class')
def confusing_resnum_mapping(test_files_sequences, test_files_structures, test_files_outputs):
    """Protein with confusing mappings between residue numbers"""
    p = Protein(ident='P08559', root_dir=test_files_outputs)
    p.load_uniprot(uniprot_id='P08559',
                   uniprot_seq_file=op.join(test_files_sequences, 'P08559.fasta'),
                   uniprot_xml_file=op.join(test_files_sequences, 'P08559.xml'),
                   download=False, set_as_representative=True)
    p.load_pdb(pdb_id='3exf_bio1', mapped_chains=['A', 'C'], pdb_file=op.join(test_files_structures, '3exf_bio1.pdb'),
               file_type='pdb', representative_chain='A', set_as_representative=True)
    p.align_seqprop_to_structprop(seqprop=p.sequences.get_by_id('P08559'),
                                  structprop=p.structures.get_by_id('3exf_bio1'),
                                  chain_id='A')
    return p

def test__map_seqprop_resnums_to_structprop_chain_index(test2, straightforward_resnum_mapping, confusing_resnum_mapping):
    assert test2._map_seqprop_resnums_to_structprop_chain_index(resnums=52, use_representatives=True) == {52: 1}

    assert straightforward_resnum_mapping._map_seqprop_resnums_to_structprop_chain_index(resnums=3,
                                                                                         use_representatives=True) == {3: 1}

    assert confusing_resnum_mapping._map_seqprop_resnums_to_structprop_chain_index(resnums=29,
                                                                                   use_representatives=True) == {29: 0}


def test_map_seqprop_resnums_to_structprop_resnums(test2, straightforward_resnum_mapping, confusing_resnum_mapping):
    assert test2.map_seqprop_resnums_to_structprop_resnums(resnums=52, use_representatives=True) == {52: 2}

    assert straightforward_resnum_mapping.map_seqprop_resnums_to_structprop_resnums(resnums=2,
                                                                                    use_representatives=True) == {2: 1}

    assert confusing_resnum_mapping.map_seqprop_resnums_to_structprop_resnums(resnums=[29,30,196],
                                                                              use_representatives=True) == {29: 0, 30: 1, 196: 167}

def test_map_structprop_resnums_to_seqprop_resnums(test2, straightforward_resnum_mapping, confusing_resnum_mapping):
    assert test2.map_structprop_resnums_to_seqprop_resnums(resnums=2, use_representatives=True) == {2: 52}

    assert straightforward_resnum_mapping.map_structprop_resnums_to_seqprop_resnums(resnums=1,
                                                                                    use_representatives=True) == {1: 2}

    assert confusing_resnum_mapping.map_structprop_resnums_to_seqprop_resnums(resnums=[0, 1, 167],
                                                                              use_representatives=True) == {0: 29, 1: 30, 167: 196}


# import os.path as op
# import unittest
#
# from ssbio.protein.structure.structprop import StructProp
#
# from ssbio.core.protein import Protein
# from ssbio.databases.kegg import KEGGProp
# from ssbio.databases.pdb import PDBProp
# from ssbio.databases.uniprot import UniProtProp
# from ssbio.protein.sequence.seqprop import SeqProp
# from ssbio.protein.structure.homology.itasser.itasserprop import ITASSERProp
#
#
# class TestProtein(unittest.TestCase):
#     """Unit tests for Protein"""
#
#     @classmethod
#     def setUpClass(self):
#         self.prot = Protein(ident='b4384')
#
#     def test_load_kegg(self):
#         new_kegg = self.prot.load_kegg('eco:b4384',
#                                        kegg_seq_file='test_files/sequences/eco-b4384.faa',
#                                        kegg_metadata_file='test_files/sequences/eco-b4384.kegg')
#         self.assertTrue(self.prot.sequences.has_id('eco:b4384'))
#         self.assertTrue(isinstance(new_kegg, KEGGProp))
#
#     def test_load_uniprot(self):
#         new_uniprot = self.prot.load_uniprot('P0ABP8',
#                                              uniprot_seq_file='test_files/sequences/P0ABP8.fasta',
#                                              uniprot_xml_file='test_files/sequences/P0ABP8.txt')
#         self.assertTrue(self.prot.sequences.has_id('P0ABP8'))
#         self.assertTrue(isinstance(new_uniprot, UniProtProp))
#
#     def test_filter_sequences(self):
#         self.prot.load_kegg('keggdummy')
#         self.prot.load_uniprot('P0OOR0')
#
#         only_kegg = self.prot.filter_sequences(KEGGProp)
#         for k in only_kegg:
#             self.assertTrue(isinstance(k, KEGGProp))
#
#         only_uniprot = self.prot.filter_sequences(UniProtProp)
#         for k in only_uniprot:
#             self.assertTrue(isinstance(k, UniProtProp))
#
#     def test_load_manual_sequence_file(self):
#         new_manual = self.prot.load_manual_sequence_file('tester', 'test_files/sequences/P0ABP8.fasta')
#         self.assertTrue(self.prot.sequences.has_id('tester'))
#         self.assertTrue(isinstance(new_manual, SeqProp))
#
#     def test_load_manual_sequence(self):
#         new_manual = self.prot.load_manual_sequence(ident='tester2', seq='ALALALAL', outdir='test_files/out/')
#         self.assertTrue(self.prot.sequences.has_id('tester2'))
#         self.assertTrue(isinstance(new_manual, SeqProp))
#
#     def test_set_representative_sequence(self):
#         self.prot.set_representative_sequence()
#         self.assertEqual(self.prot.representative_sequence.id, 'P0ABP8')
#         self.assertTrue(self.prot.representative_sequence.num_pdbs > 0)
#
#     def test_align_sequences_to_representative(self):
#         # TODO: unittests for this method
#         # this is only true when you dont add more sequences
#         # self.assertEqual(len(self.prot.representative_sequence.alignments), len(self.prot.sequences))
#         pass
#
#     def test_load_pdb(self):
#         new_uniprot = self.prot.load_uniprot('P0ABP8',
#                                              uniprot_seq_file='test_files/sequences/P0ABP8.fasta',
#                                              uniprot_xml_file='test_files/sequences/P0ABP8.txt',
#                                              set_as_representative=True)
#
#         newpdb = self.prot.load_pdb('1ecp', pdb_file='test_files/structures/1ecp.pdb', file_type='pdb')
#         self.assertTrue(self.prot.structures.has_id('1ecp'))
#         self.assertTrue(isinstance(newpdb, PDBProp))
#
#     def test_load_itasser_folder(self):
#         self.assertRaises(OSError, self.prot.load_itasser_folder, 'sdfsdfsdf', itasser_folder='sdfsdfsdf/structures/P9WG73')
#
#         newitasser = self.prot.load_itasser_folder('P9WG73', itasser_folder='test_files/structures/P9WG73')
#         self.assertTrue(self.prot.structures.has_id('P9WG73'))
#         self.assertTrue(isinstance(newitasser, ITASSERProp))
#         self.assertEqual(newitasser.structure_file, 'model1.pdb')
#
#         # newitasser.copy_results(copy_to_dir='test_files/structures/test_out/', rename_model_to='haha', force_rerun=True)
#
#     def test_load_homology_model(self):
#         newprot = self.prot.load_pdb_file('DEOD', pdb_file='test_files/structures/DEOD_ECOLI_model1.pdb')
#         self.assertTrue(self.prot.structures.has_id('DEOD'))
#         self.assertTrue(isinstance(newprot, StructProp))
#
#     def test_num_structures(self):
#         self.assertNotEqual(self.prot.num_structures, 0)
#
#     def test_num_structures_experimental(self):
#         self.assertNotEqual(self.prot.num_structures_experimental, 0)
#
#     def test_num_structures_homology(self):
#         self.assertNotEqual(self.prot.num_structures_homology, 0)
#
#     def test_set_representative_structure(self):
#         prot_to_set = self.prot.structures.get_by_id('1ecp')
#         self.prot.set_representative_structure(seq_outdir=op.join('test_files','out'),
#                                                struct_outdir=op.join('test_files','out'),
#                                                pdb_file_type='cif')
#         self.assertEqual('1ecp-A', self.prot.representative_structure.id)