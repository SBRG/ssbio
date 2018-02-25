import os.path as op
import pytest
from ssbio.pipeline.gempro import GEMPRO
from ssbio.core.genepro import GenePro
from ssbio.databases.kegg import KEGGProp
from ssbio.databases.uniprot import UniProtProp
from ssbio.databases.pdb import PDBProp
from cobra.core import Model, DictList


@pytest.fixture(scope='class')
def gempro_empty():
    """GEMPRO with ID"""
    return GEMPRO(gem_name='test_id')


@pytest.fixture(scope='class')
def gempro_with_dir(test_files_tempdir):
    """GEMPRO with ID and set root directory"""
    return GEMPRO(gem_name='test_id_dir', root_dir=test_files_tempdir)


@pytest.fixture(scope='class')
def gempro_with_dir_mini_json(test_files_tempdir, test_gem_small_json):
    """GEMPRO with ID + GEM (mini JSON)"""
    return GEMPRO(gem_name='test_id_dir_json', root_dir=test_files_tempdir,
                  gem_file_path=test_gem_small_json, gem_file_type='json')


@pytest.fixture(scope='class')
def gempro_with_dir_genes(test_files_tempdir, list_of_genes_ecoli):
    """GEMPRO with ID + list of gene IDs"""
    return GEMPRO(gem_name='test_id_geneids', root_dir=test_files_tempdir, genes_list=list_of_genes_ecoli)


@pytest.fixture(scope='class')
def gempro_with_dir_genes_seqs(test_files_tempdir, dict_of_genes_seqs_ecoli):
    """GEMPRO with ID + dictionary of gene IDs + protein sequences"""
    return GEMPRO(gem_name='test_id_geneseqs', root_dir=test_files_tempdir, genes_and_sequences=dict_of_genes_seqs_ecoli)


@pytest.fixture(scope='class')
def gempro_with_dir_genes_fasta(test_files_tempdir, test_fasta_file_multseqs):
    """GEMPRO with ID + FASTA file"""
    return GEMPRO(gem_name='test_id_fasta', root_dir=test_files_tempdir,
                  genome_path=test_fasta_file_multseqs)


def test_directory_creation(gempro_with_dir, test_files_tempdir):
    check_dir = op.join(test_files_tempdir, gempro_with_dir.id)
    assert op.exists(check_dir)


class TestGemproWithDirMiniJson():
    """Tests for the gempro_with_dir_mini_json fixture"""

    @pytest.mark.run(order=1)
    def test_init(self, gempro_with_dir_mini_json):
        """Test initialization of a GEMPRO"""
        assert op.exists(gempro_with_dir_mini_json.base_dir)
        assert op.exists(gempro_with_dir_mini_json.data_dir)
        assert op.exists(gempro_with_dir_mini_json.genes_dir)

        # Test that we can access the model using .model
        assert isinstance(gempro_with_dir_mini_json.model, Model)

    @pytest.mark.run(order=2)
    def test_gene_info(self, gempro_with_dir_mini_json):
        """Test that genes is a DictList, contains GenePro objects, and directories were made"""
        assert isinstance(gempro_with_dir_mini_json.genes, DictList)
        for g in gempro_with_dir_mini_json.genes:
            assert isinstance(g, GenePro)
            assert op.exists(g.gene_dir)
            assert op.exists(g.protein.protein_dir)

    @pytest.mark.run(order=3)
    def test_add_remove_gene_ids(self, gempro_with_dir_mini_json):
        """Test adding and removing genes to a GEMPRO, mainly to double check the genes convenience attribute"""
        add_these = ['b0002', 'b2092']
        gempro_with_dir_mini_json.add_gene_ids(add_these)
        for x in add_these:
            assert gempro_with_dir_mini_json.genes.has_id(x)
            assert isinstance(gempro_with_dir_mini_json.genes.get_by_id(x), GenePro)

        remove_these = ['b4395', 'b3612', 'b2417', 'b2779', 'b2925', 'b2926', 'b3916', 'b3919', 'b4025', 'b2097']
        for x in remove_these:
            gempro_with_dir_mini_json.genes.remove(gempro_with_dir_mini_json.genes.get_by_id(x))
            assert not gempro_with_dir_mini_json.genes.has_id(x)

        assert len(gempro_with_dir_mini_json.genes) == len(gempro_with_dir_mini_json.model.genes)

    @pytest.mark.run(order=4)
    def test_kegg_mapping_and_metadata(self, gempro_with_dir_mini_json):
        """Test KEGG mapping"""

        gempro_with_dir_mini_json.kegg_mapping_and_metadata(kegg_organism_code='eco')

        for g in gempro_with_dir_mini_json.genes:

            kegg_mappings = g.protein.filter_sequences(KEGGProp)

            # Test that KEGGProp objects have not been added to the sequences if missing
            if g.id in gempro_with_dir_mini_json.missing_kegg_mapping:
                assert len(kegg_mappings) == 0
                continue

            # Test that KEGGProp objects have been added to the sequences
            assert len(kegg_mappings) > 0

            # Test that KEGG sequence files have been saved
            for k in kegg_mappings:
                assert op.exists(k.sequence_path)

    @pytest.mark.run(order=5)
    def test_uniprot_mapping_and_metadata(self, gempro_with_dir_mini_json):
        """Test UniProt mapping"""

        gempro_with_dir_mini_json.uniprot_mapping_and_metadata(model_gene_source='ENSEMBLGENOME_ID')

        for g in gempro_with_dir_mini_json.genes:

            uniprot_mappings = g.protein.filter_sequences(UniProtProp)

            # Test that UniProtProp objects have not been added to the sequences if missing
            if g.id in gempro_with_dir_mini_json.missing_uniprot_mapping:
                assert len(uniprot_mappings) == 0
                continue

            # Test that UniProtProp objects have been added to the sequences
            assert len(uniprot_mappings) > 0

            # Test that UniProt XML files have been saved
            for u in uniprot_mappings:
                assert op.exists(u.metadata_path)

    @pytest.mark.run(order=6)
    def test_set_representative_sequence(self, gempro_with_dir_mini_json):
        """Test representative sequence setting"""

        gempro_with_dir_mini_json.set_representative_sequence()

        # Test proteins that shouldn't map don't have a representative sequence
        will_not_map = ['b2092']
        for g_id in will_not_map:
            assert g_id in gempro_with_dir_mini_json.missing_representative_sequence
            assert not gempro_with_dir_mini_json.genes_with_a_representative_sequence.has_id(g_id)

    @pytest.mark.run(order=7)
    def test_map_uniprot_to_pdb(self, gempro_with_dir_mini_json):
        """Test mapping of UniProt IDs to PDB IDs"""

        # Mapped to PDB successfully as of 2018-02-24
        mapped_to_pdb = ['b0755', 'b0875', 'b1101', 'b1676', 'b1723', 'b1779', 'b1817', 'b2133', 'b2415', 'b2416']

        gempro_with_dir_mini_json.map_uniprot_to_pdb(seq_ident_cutoff=0.0)

        for g_id in mapped_to_pdb:
            g = gempro_with_dir_mini_json.genes.get_by_id(g_id)
            assert g.protein.num_structures_experimental > 0

            for s in g.protein.get_experimental_structures():
                # Test that PDBProp objects have been added to protein structures list
                assert isinstance(s, PDBProp)
                # Test that their structures have NOT been downloaded yet
                assert not s.structure_file

    @pytest.mark.run(order=8)
    def test_blast_seqs_to_pdb(self, gempro_with_dir_mini_json):
        """Test BLASTing sequences to PDB"""

        # BLASTed to PDB successfully as of 2018-02-24
        additional_blasted_to_pdb = ['b1380']

        gempro_with_dir_mini_json.blast_seqs_to_pdb(seq_ident_cutoff=0.5, evalue=0.0001, all_genes=False)

        for g_id in additional_blasted_to_pdb:
            p = gempro_with_dir_mini_json.genes.get_by_id(g_id).protein
            assert p.num_structures_experimental > 0

            for s in p.get_experimental_structures():
                # Test that PDBProp objects have been added to protein structures list
                assert isinstance(s, PDBProp)
                # Test that their structures have NOT been downloaded yet
                assert not s.structure_file

    @pytest.mark.run(order=9)
    def test_set_representative_structure(self, gempro_with_dir_mini_json):
        """Test setting the representative structure"""

        gempro_with_dir_mini_json.set_representative_structure(engine='biopython', always_use_homology=False,
                                                               rez_cutoff=0.0, seq_ident_cutoff=0.5,
                                                               allow_missing_on_termini=0.2,
                                                               allow_mutants=True, allow_deletions=False,
                                                               allow_insertions=False, allow_unresolved=True,
                                                               skip_large_structures=False, clean=True)

        # List of genes with representative structures as of 2018-02-24
        genes_with_repstructs = ['b0755', 'b0875', 'b1723', 'b1779', 'b2415', 'b2416']

        # TODO: should we test the large structure setting here?
        genes_with_large_struct = []

        for g_id in genes_with_repstructs:
            p = gempro_with_dir_mini_json.genes.get_by_id(g_id).protein
            assert p.representative_structure
            assert p.structures.has_id(p.representative_structure.id)
            assert op.exists(p.representative_structure.structure_path)