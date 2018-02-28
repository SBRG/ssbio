import os.path as op
import pytest
import ssbio.io
import shutil
import pandas as pd
from ssbio.pipeline.atlas2 import ATLAS2


@pytest.fixture(scope='class')
def reference_gempro(test_files_gempro):
    gp = ssbio.io.load_json(op.join(test_files_gempro, 'test_id_dir_json', 'model', 'mini_gp.json'))
    gp.root_dir = test_files_gempro
    return gp


@pytest.fixture(scope='class')
def test_strains():
    return ['585395.4', '1169323.3', '1068619.3', '1005537.3', '868163.3']


@pytest.fixture(scope='class')
def strains_to_fasta(test_strains, test_files_atlas):
    d = {}
    for s in test_strains:
        d[s] = op.join(test_files_atlas, '{}.PATRIC.faa'.format(s))
    return d


@pytest.fixture(scope='class')
def premade_orthology(test_files_atlas):
    """Custom altered orthology matrix"""
    return op.join(test_files_atlas, 'mini_orth_matrix.csv')


@pytest.fixture(scope='class')
def my_atlas(reference_gempro, test_files_tempdir):
    return ATLAS2(atlas_name='test_atlas', root_dir=test_files_tempdir,
                  reference_gempro=reference_gempro, reference_genome_path=reference_gempro.genome_path)


class TestATLAS():
    @pytest.mark.run(order=1)
    def test_init(self, my_atlas):
        assert op.exists(my_atlas.base_dir)
        assert op.exists(my_atlas.model_dir)
        assert op.exists(my_atlas.data_dir)
        assert op.exists(my_atlas.sequences_dir)
        assert op.exists(my_atlas.sequences_by_gene_dir)
        assert op.exists(my_atlas.sequences_by_organism_dir)

    @pytest.mark.run(order=2)
    def test_load_strains(self, my_atlas, strains_to_fasta):

        my_atlas.load_strains(strains_to_fasta)

        assert len(my_atlas.strains) == len(strains_to_fasta)
        for s in my_atlas.strains:
            assert s.id in strains_to_fasta
            # Genome paths should be linked
            assert s.genome_path == strains_to_fasta[s.id]
            # No genes stored yet
            assert len(s.genes) == 0

    @pytest.mark.run(order=3)
    def test_get_orthology_matrix(self, my_atlas, premade_orthology):
        # Have to copy the premade DF to the atlas data_dir, just check if it was copied
        shutil.copy(premade_orthology, my_atlas.data_dir)
        assert op.exists(op.join(my_atlas.data_dir, op.basename(premade_orthology)))

        my_atlas.get_orthology_matrix(outfile=op.basename(premade_orthology), force_rerun=False)
        assert isinstance(my_atlas.df_orthology_matrix, pd.DataFrame)

    @pytest.mark.run(order=4)
    def test_filter_genes_and_strains(self, my_atlas, test_strains):

        for g in my_atlas.reference_gempro.genes:
            assert g.functional

        functional_genes = ['b0755', 'b0875', 'b1101', 'b1380', 'b1621', 'b1676', 'b1723', 'b1773', 'b1779', 'b1817',
                            'b1818', 'b1819', 'b1854', 'b2097', 'b2133', 'b2415', 'b2416', 'b2975', 'b2987', 'b3493',
                            'b3603', 'b2092', 'b0002']

        my_atlas.filter_genes_and_strains(custom_keep_strains=test_strains[:4])

        assert len(my_atlas.strains) == 4

        for g in my_atlas.reference_gempro.genes:
            if g.id in functional_genes:
                assert g.functional
            else:
                assert not g.functional

    @pytest.mark.run(order=5)
    def test_build_strain_specific_models(self, my_atlas):

        my_atlas.build_strain_specific_models(force_rerun=True)

        for s in my_atlas.strains:
            assert op.exists(op.join(my_atlas.model_dir, '{}_gp.pckl'.format(s.id)))
            assert len(s.genes) == len(my_atlas.reference_gempro.genes)
            assert len(s.functional_genes) < len(s.genes)

            for g in my_atlas.reference_gempro.functional_genes:
                if s.id == '585395.4' and g.id == 'b1773':
                    assert not s.genes.get_by_id(g.id).functional
                elif s.id == '1169323.3' and g.id == 'b2133':
                    assert not s.genes.get_by_id(g.id).functional
                elif s.id == '1068619.3' and g.id == 'b1779':
                    assert not s.genes.get_by_id(g.id).functional
                elif s.id == '1005537.3' and (g.id == 'b2975' or g.id == 'b2987'):
                    assert not s.genes.get_by_id(g.id).functional
                else:
                    assert s.genes.get_by_id(g.id).functional

    @pytest.mark.run(order=6)
    def test_load_sequences_to_strains(self, my_atlas):

        my_atlas.load_sequences_to_strains(force_rerun=True)

        for s in my_atlas.strains:
            assert op.exists(op.join(my_atlas.model_dir, '{}_gp_withseqs.pckl'.format(s.id)))
            for g in s.functional_genes:
                assert g.protein.representative_sequence

    @pytest.mark.run(order=7)
    def test_load_sequences_to_reference(self, my_atlas):

        my_atlas.load_sequences_to_reference(force_rerun=True)

        for g in my_atlas.reference_gempro.functional_genes:
            assert op.exists(op.join(my_atlas.sequences_by_gene_dir, '{}_protein_withseqs.pckl'.format(g.id)))
            for s in my_atlas.strains:
                check_id = '{}_{}'.format(s.genes.get_by_id(g.id).id, s.id)
                if s.genes.get_by_id(g.id).functional:
                    assert g.protein.sequences.has_id(check_id)
                else:
                    assert not g.protein.sequences.has_id(check_id)

    @pytest.mark.run(order=8)
    def test_align_orthologous_genes_pairwise(self, my_atlas):

        my_atlas.align_orthologous_genes_pairwise(engine='biopython', force_rerun=True)

        for g in my_atlas.reference_gempro.functional_genes:
            if g.protein.representative_sequence:
                assert op.exists(op.join(my_atlas.sequences_by_gene_dir, '{}_protein_withseqs_aln.pckl'.format(g.id)))
                for s in my_atlas.strains:
                    check_id = '{}_{}_{}'.format(g.id, s.genes.get_by_id(g.id).id, s.id)
                    if s.genes.get_by_id(g.id).functional:
                        assert g.protein.sequence_alignments.has_id(check_id)
                    else:
                        assert not g.protein.sequence_alignments.has_id(check_id)
            else:
                assert len(g.protein.sequence_alignments) == 0