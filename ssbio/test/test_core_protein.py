import pytest
from ssbio.core.protein import Protein
import os.path as op
from pathlib import Path

from ssbio.databases.kegg import KEGGProp
from ssbio.databases.pdb import PDBProp
from ssbio.databases.uniprot import UniProtProp
from ssbio.protein.sequence.seqprop import SeqProp
from ssbio.protein.structure.homology.itasser.itasserprop import ITASSERProp


@pytest.fixture(scope='class')
def FixtProtein(test_files_sequences, test_files_structures, test_files_outputs):
    """General Protein"""
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
def FixtProtein_straightforwardmapping(test_files_sequences, test_files_structures, test_files_outputs):
    """Protein with straightforward mappings between residue numbers"""
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
def FixtProtein_confusingmapping(test_files_sequences, test_files_structures, test_files_outputs):
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


def test_map_seqprop_resnums_to_structprop_chain_index(FixtProtein, FixtProtein_straightforwardmapping, FixtProtein_confusingmapping):
    assert FixtProtein._map_seqprop_resnums_to_structprop_chain_index(
        resnums=52,
        use_representatives=True) == {52: 1}

    assert FixtProtein_straightforwardmapping._map_seqprop_resnums_to_structprop_chain_index(
        resnums=3,
        use_representatives=True) == {3: 1}

    assert FixtProtein_confusingmapping._map_seqprop_resnums_to_structprop_chain_index(
        resnums=29,
        use_representatives=True) == {29: 0}


def test_map_seqprop_resnums_to_structprop_resnums(FixtProtein, FixtProtein_straightforwardmapping, FixtProtein_confusingmapping):
    assert FixtProtein.map_seqprop_resnums_to_structprop_resnums(
        resnums=52,
        use_representatives=True) == {52: 2}

    assert FixtProtein_straightforwardmapping.map_seqprop_resnums_to_structprop_resnums(
        resnums=2,
        use_representatives=True) == {2: 1}

    assert FixtProtein_confusingmapping.map_seqprop_resnums_to_structprop_resnums(
        resnums=[29,30,196],
        use_representatives=True) == {29: 0, 30: 1, 196: 167}


def test_map_structprop_resnums_to_seqprop_resnums(FixtProtein, FixtProtein_straightforwardmapping, FixtProtein_confusingmapping):
    assert FixtProtein.map_structprop_resnums_to_seqprop_resnums(
        resnums=2,
        use_representatives=True) == {2: 52}

    assert FixtProtein_straightforwardmapping.map_structprop_resnums_to_seqprop_resnums(
        resnums=1,
        use_representatives=True) == {1: 2}

    assert FixtProtein_confusingmapping.map_structprop_resnums_to_seqprop_resnums(
        resnums=[0, 1, 167],
        use_representatives=True) == {0: 29, 1: 30, 167: 196}


@pytest.fixture(scope='class')
def FixtProtein_comprehensive(test_files_sequences, test_files_structures, test_files_outputs):
    """General Protein 2"""
    p = Protein(ident='b4384', root_dir=test_files_outputs)

    kegg_seq_file = Path(test_files_sequences) / "eco-b4384.faa"
    kegg_metadata_file = Path(test_files_sequences) / "eco-b4384.kegg"
    new_kegg = p.load_kegg(
        'eco:b4384',
        kegg_seq_file=kegg_seq_file,
        kegg_metadata_file=kegg_metadata_file)

    uniprot_seq_file = Path(test_files_sequences) / 'P0ABP8.fasta'
    uniprot_xml_file = Path(test_files_sequences) / 'P0ABP8.xml'
    new_uniprot = p.load_uniprot(
        'P0ABP8',
        uniprot_seq_file=uniprot_seq_file,
        uniprot_xml_file=uniprot_xml_file)

    pdb_file = Path(test_files_structures) / '1ecp.pdb'
    newpdb = p.load_pdb(
        '1ecp',
        pdb_file=pdb_file,
        file_type='pdb')

    return p, new_kegg, new_uniprot, newpdb


def test_load_kegg(FixtProtein_comprehensive):
    assert FixtProtein_comprehensive[0].sequences.has_id('eco:b4384')
    assert isinstance(FixtProtein_comprehensive[1], KEGGProp)


def test_load_uniprot(FixtProtein_comprehensive):
    assert FixtProtein_comprehensive[0].sequences.has_id('P0ABP8')
    assert isinstance(FixtProtein_comprehensive[2], UniProtProp)


def test_filter_sequences(FixtProtein_comprehensive):
    FixtProtein_comprehensive[0].load_kegg('keggdummy')
    FixtProtein_comprehensive[0].load_uniprot('P0OOR0')

    only_kegg = FixtProtein_comprehensive[0].filter_sequences(KEGGProp)
    for k in only_kegg:
        assert isinstance(k, KEGGProp)

    only_uniprot = FixtProtein_comprehensive[0].filter_sequences(UniProtProp)
    for k in only_uniprot:
        assert isinstance(k, UniProtProp)


def test_load_manual_sequence_file(FixtProtein_comprehensive, test_files_sequences):
    new_manual = FixtProtein_comprehensive[0].load_manual_sequence_file('tester', Path(test_files_sequences) / 'P0ABP8.fasta')
    assert FixtProtein_comprehensive[0].sequences.has_id('tester')
    assert isinstance(new_manual, SeqProp)


def test_load_manual_sequence(FixtProtein_comprehensive, test_files_outputs):
    new_manual = FixtProtein_comprehensive[0].load_manual_sequence(ident='tester2', seq='ALALALAL', outdir=test_files_outputs)
    assert FixtProtein_comprehensive[0].sequences.has_id('tester2')
    assert isinstance(new_manual, SeqProp)


def test_set_representative_sequence(FixtProtein_comprehensive):
    FixtProtein_comprehensive[0].set_representative_sequence()
    assert FixtProtein_comprehensive[0].representative_sequence.id == 'P0ABP8'
    assert FixtProtein_comprehensive[0].representative_sequence.num_pdbs > 0


def test_align_sequences_to_representative(FixtProtein_comprehensive):
    # TODO: unittests for this method
    # this is only true when you dont add more sequences
    # self.assertEqual(len(FixtProtein_comprehensive.representative_sequence.alignments), len(FixtProtein_comprehensive.sequences))
    pass


def test_load_pdb(FixtProtein_comprehensive):
    assert FixtProtein_comprehensive[0].structures.has_id('1ecp')
    assert isinstance(FixtProtein_comprehensive[3], PDBProp)

#
# def test_load_itasser_folder(FixtProtein_comprehensive):
#     self.assertRaises(OSError, FixtProtein_comprehensive.load_itasser_folder, 'sdfsdfsdf', itasser_folder='sdfsdfsdf/structures/P9WG73')
#
#     newitasser = FixtProtein_comprehensive.load_itasser_folder('P9WG73', itasser_folder='test_files/structures/P9WG73')
#     assert FixtProtein_comprehensive.structures.has_id('P9WG73')
#     assert isinstance(newitasser, ITASSERProp)
#     self.assertEqual(newitasser.structure_file, 'model1.pdb')
#
#     # newitasser.copy_results(copy_to_dir='test_files/structures/test_out/', rename_model_to='haha', force_rerun=True)

# def test_load_homology_model(FixtProtein_comprehensive):
#     newprot = FixtProtein_comprehensive.load_pdb_file('DEOD', pdb_file='test_files/structures/DEOD_ECOLI_model1.pdb')
#     assert FixtProtein_comprehensive.structures.has_id('DEOD')
#     assert isinstance(newprot, StructProp)
#
# def test_num_structures(FixtProtein_comprehensive):
#     self.assertNotEqual(FixtProtein_comprehensive.num_structures, 0)
#
# def test_num_structures_experimental(FixtProtein_comprehensive):
#     self.assertNotEqual(FixtProtein_comprehensive.num_structures_experimental, 0)
#
# def test_num_structures_homology(FixtProtein_comprehensive):
#     self.assertNotEqual(FixtProtein_comprehensive.num_structures_homology, 0)
#
# def test_set_representative_structure(FixtProtein_comprehensive):
#     prot_to_set = FixtProtein_comprehensive.structures.get_by_id('1ecp')
#     FixtProtein_comprehensive.set_representative_structure(seq_outdir=op.join('test_files','out'),
#                                            struct_outdir=op.join('test_files','out'),
#                                            pdb_file_type='cif')
#     self.assertEqual('1ecp-A', FixtProtein_comprehensive.representative_structure.id)