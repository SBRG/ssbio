import unittest
import tempfile
import os.path as op
from ssbio.gempro.pipeline import GEMPRO
import cobra
from ssbio.utils import Date
date = Date()

class TestGEMPRO(unittest.TestCase):
    """Unit tests for the GEMPRO pipeline
    """

    @classmethod
    def setUpClass(self):
        self.GEM_NAME = 'ecoli_test'
        self.tempdir = tempfile.TemporaryDirectory()
        self.ROOT_DIR = self.tempdir.name
        self.my_gempro = GEMPRO(self.GEM_NAME, self.ROOT_DIR)
        self.my_gempro.prep_folders()

        gem_file = 'test_files/Ec_core_flux1.xml'
        self.my_gempro.load_model(gem_file=gem_file, file_type='sbml')

        # remove some genes from this model so testing doesn't take too long
        remove_these = ['b0118', 'b4025', 'b4153', 'b2925', 'b3919', 'b3738', 'b3732', 'b0726', 'b1602', 'b1101',
                        'b3236', 'b0728', 'b1603', 'b2926', 'b0432', 'b3735', 'b0474', 'b4090', 'b3731', 'b0767',
                        'b3737', 'b0724', 'b0008', 'b3403', 'b1779', 'b0727', 'b3956', 'b1676', 'b2935', 'b4232',
                        'b0430', 'b2914', 'b0723', 'b0722', 'b0729', 'b0720', 'b0116', 'b4152', 'b2415', 'b0429',
                        'b0621', 'b0904', 'b1241', 'b1276', 'b1380', 'b1416', 'b1479', 'b1612', 'b1702', 'b2417',
                        'b1852', 'b1854', 'b2276', 'b2277', 'b2278', 'b2279', 'b2280', 'b2281', 'b2282', 'b2283',
                        'b2284', 'b2285', 'b2286', 'b2287', 'b2288', 'b2297', 'b2463', 'b2465', 'b2587', 'b2975',
                        'b3117', 'b3386', 'b3528', 'b3603', 'b3739', 'b3951', 'b3952', 'b4122', 'b4395', 'b0721',
                        'b2029', 'b1136', 'b4015', 'b4014', 'b2976', 'b0114', 'b0115', 'b3916', 'b1723', 'b0755'
                        ]
        for x in remove_these:
            self.my_gempro.genes.remove(x)

        # test the manual adding of any gene ID
        add_these = ['b0002', 'b0003', 'b0004', 'b2092']
        self.my_gempro.genes.extend(add_these)

        self.will_not_map = ['b1417', 'b2092']

    def test_gene_content(self):
        # test if there are only these genes remaining
        remaining = ['b2296', 'b3734', 'b3733', 'b3736', 'b0431', 'b2779', 'b4151', 'b4154', 'b1611',  'b2416',
                     'b0002', 'b0003', 'b0004', 'b1417', 'b2092']
        self.assertCountEqual(self.my_gempro.genes, remaining)

    def test_prep_folders(self):
        # test if folders were made
        self.assertTrue(op.exists(op.join(self.ROOT_DIR, self.GEM_NAME)))
        folders = ['data', 'notebooks', 'figures', 'structure_files', 'structure_files/by_gene', 'sequence_files']
        for f in folders:
            self.assertTrue(op.exists(op.join(self.ROOT_DIR, self.GEM_NAME, f)))

    def test_load_model(self):
        # test that we can access the model using .model
        self.assertIsInstance(self.my_gempro.model, cobra.Model)

        # test original model attributes
        self.assertEqual(len(self.my_gempro.model.reactions), 77)
        self.assertEqual(len(self.my_gempro.model.metabolites), 63)
        self.assertEqual(len(self.my_gempro.model.genes), 101)

    def test_kegg_mapping_and_metadata(self):
        kegg_organism_code = 'eco'
        outfile_df = self.my_gempro.kegg_mapping_and_metadata(kegg_organism_code=kegg_organism_code, force_rerun=True)

        # test if mapping dataframe was created
        self.assertTrue(op.exists(op.join(self.my_gempro.data, outfile_df)))

        # test if sequences and metadata were downloaded
        for g in self.my_gempro.genes:
            gene_structure_folder = op.join(self.my_gempro.seq_files, g)

            # if the gene did not map to a KEGG ID, a folder will still be made for it
            # it will still be in the KEGG DF, just with empty entries
            # it will not have anything downloaded for it
            self.assertTrue(op.exists(gene_structure_folder))
            self.assertTrue(g in self.my_gempro.kegg_df.m_gene.tolist())
            if g in self.will_not_map:
                self.assertTrue(g in self.my_gempro.kegg_missing)
                self.assertFalse(op.exists(op.join(gene_structure_folder, '{}-{}.faa'.format(kegg_organism_code, g))))
            else:
                self.assertTrue(op.exists(op.join(gene_structure_folder, '{}-{}.faa'.format(kegg_organism_code, g))))
                self.assertTrue(op.exists(op.join(gene_structure_folder, '{}-{}.kegg'.format(kegg_organism_code, g))))

    def test_uniprot_mapping_and_metadata(self):
        outfile_df = self.my_gempro.uniprot_mapping_and_metadata('ENSEMBLGENOME_ID', force_rerun=True)

        # test if mapping dataframe was created
        self.assertTrue(op.exists(op.join(self.my_gempro.data, outfile_df)))

        # test if sequences and metadata were downloaded
        for g in self.my_gempro.genes:
            gene_structure_folder = op.join(self.my_gempro.seq_files, g)

            # if the gene did not map to a UniProt ID, a folder will still be made for it
            # it will still be in the UniProt DF, just with empty entries
            # it will not have anything downloaded for it
            self.assertTrue(op.exists(gene_structure_folder))
            self.assertTrue(g in self.my_gempro.uniprot_df.m_gene.tolist())
            if g in self.will_not_map:
                self.assertTrue(g in self.my_gempro.uniprot_missing)
                self.assertFalse(op.exists(op.join(gene_structure_folder, '{}.fasta'.format(u))))
                self.assertFalse(op.exists(op.join(gene_structure_folder, '{}.txt'.format(u))))
            else:
                uniprots = self.my_gempro.uniprot_df[self.my_gempro.uniprot_df.m_gene == g].u_uniprot_acc.unique().tolist()
                for u in uniprots:
                    self.assertTrue(op.exists(op.join(gene_structure_folder, '{}.fasta'.format(u))))
                    self.assertTrue(op.exists(op.join(gene_structure_folder, '{}.txt'.format(u))))

    def test_consolidate_mappings(self):
        pass

    @classmethod
    def tearDownClass(self):
        self.tempdir.cleanup()