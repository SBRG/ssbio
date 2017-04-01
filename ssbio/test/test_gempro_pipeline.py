import shutil
import os.path as op
import unittest
import cobra
import pandas as pd
import six
from ssbio.pipeline.gempro import GEMPRO
from ssbio.core.genepro import GenePro
from ssbio.utils import Date
from ssbio.databases.kegg import KEGGProp

date = Date()

class TestGEMPRO(unittest.TestCase):
    """Unit tests for the GEMPRO pipeline
    """

    @classmethod
    def setUpClass(self):
        self.GEM_NAME = 'ecoli_test'
        self.ROOT_DIR = op.join('test_files', 'out')
        gem_file = op.join('test_files', 'Ec_core_flux1.xml')

        self.my_gempro = GEMPRO(gem_name=self.GEM_NAME, root_dir=self.ROOT_DIR, gem_file_path=gem_file, gem_file_type='sbml')
        self.base_dir = self.my_gempro.base_dir

        # Remove some genes from this model so testing doesn't take too long
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
            self.my_gempro.genes.remove(self.my_gempro.genes.get_by_id(x))

        # Test the manual adding of any gene ID
        add_these = ['b0002', 'b0003', 'b0004', 'b2092']
        self.my_gempro.add_genes_by_id(add_these)

        self.kegg_will_not_map = ['b1417', 'b2092']
        self.uniprot_will_not_map = ['b1417', 'b2092']

        print(self.base_dir)

    def test_add_genes_by_id(self):
        add_these = ['b0002', 'b0003', 'b0004', 'b2092']
        for x in add_these:
            self.assertTrue(self.my_gempro.genes.has_id(x))

    def test_gene_content(self):
        # Test that genes were correctly removed from self.genes
        remaining = ['b2296', 'b3734', 'b3733', 'b3736', 'b0431', 'b2779', 'b4151', 'b4154', 'b1611',  'b2416',
                     'b0002', 'b0003', 'b0004', 'b1417', 'b2092']
        six.assertCountEqual(self, [x.id for x in self.my_gempro.genes], remaining)

    def test_prep_folders(self):
        # Test if folders were made
        self.assertTrue(op.exists(op.join(self.ROOT_DIR, self.GEM_NAME)))
        folders = ['data', 'structures', 'sequences']
        for f in folders:
            self.assertTrue(op.exists(op.join(self.ROOT_DIR, self.GEM_NAME, f)))

    def test_load_model(self):
        # Test that we can access the model using .model
        self.assertIsInstance(self.my_gempro.model, cobra.Model)

        # Test original model attributes
        self.assertEqual(len(self.my_gempro.model.reactions), 77)
        self.assertEqual(len(self.my_gempro.model.metabolites), 63)
        # genes should be removed from model
        self.assertEqual(len(self.my_gempro.genes), 15)

    def test_gene_info(self):
        # Test that self.genes is a DictList and contains GenePro objects
        self.assertIsInstance(self.my_gempro.genes, cobra.DictList)
        for g in self.my_gempro.genes:
            self.assertIsInstance(g, GenePro)

    def test_kegg_mapping_and_metadata(self):
        """Test that KEGG mapping did these things:
        1. Create the df_kegg_metadata attribute
        2. Save individual sequence and metadata files per gene in the sequence_files folder
        3. Save KEGG related info into each GeneInfo object
        """
        kegg_organism_code = 'eco'
        self.my_gempro.kegg_mapping_and_metadata(kegg_organism_code=kegg_organism_code)
        # print(self.my_gempro.df_kegg_metadata)

        # 1.
        # Test if attribute is a DataFrame and the number of entries in it is equal to the number of genes
        self.assertIsInstance(self.my_gempro.df_kegg_metadata, pd.DataFrame)
        self.assertEqual(len(self.my_gempro.df_kegg_metadata.gene.unique()), len(self.my_gempro.genes))

        ## Tests per gene
        for gene in self.my_gempro.genes:
            g = gene.id

            # 2.
            # Test that gene metadata was downloaded
            # If the gene did not map to a KEGG ID, a folder will still be made for it
            # It will still be in the KEGG DF, just with empty entries
            # It will not have anything downloaded for it
            gene_structure_folder = op.join(self.my_gempro.sequence_dir, g)
            self.assertTrue(op.exists(gene_structure_folder))
            self.assertTrue(g in self.my_gempro.df_kegg_metadata.gene.tolist())

            if g in self.kegg_will_not_map:
                self.assertTrue(g in self.my_gempro.missing_kegg_mapping)
                self.assertFalse(op.exists(op.join(gene_structure_folder, '{}-{}.faa'.format(kegg_organism_code, g))))
            else:
                self.assertTrue(op.exists(op.join(gene_structure_folder, '{}-{}.faa'.format(kegg_organism_code, g))))
                self.assertTrue(op.exists(op.join(gene_structure_folder, '{}-{}.kegg'.format(kegg_organism_code, g))))
                # 3.
                # Test that the sequence_properties attribute has KEGG information in it
                keggs = gene.protein.filter_sequences(KEGGProp)
                for k in keggs:
                    self.assertTrue(k.seq_len > 0)

    def test_uniprot_mapping_and_metadata(self):
        """Test that UniProt mapping did these things:
        1. Create the df_uniprot_metadata attribute
        2. Save individual sequence and metadata files per gene in the sequence_files folder
        3. Save UniProt related info into each GeneInfo object
        """
        self.my_gempro.uniprot_mapping_and_metadata(model_gene_source='ENSEMBLGENOME_ID')
        # print(self.my_gempro.df_uniprot_metadata)

        # 1.
        # Test if attribute is a DataFrame and the number of entries in it is equal to the number of genes
        self.assertIsInstance(self.my_gempro.df_uniprot_metadata, pd.DataFrame)
        self.assertEqual(len(self.my_gempro.df_uniprot_metadata.gene.unique()), len(self.my_gempro.genes))

        ## Tests per gene
        for gene in self.my_gempro.genes:
            g = gene.id

            # 2.
            # Test that gene metadata was downloaded
            # If the gene did not map to a UniProt ID, a folder will still be made for it
            # It will still be in the UniProt DF, just with empty entries
            # It will not have anything downloaded for it
            gene_structure_folder = op.join(self.my_gempro.sequence_dir, g)
            self.assertTrue(op.exists(gene_structure_folder))
            self.assertTrue(g in self.my_gempro.df_uniprot_metadata.gene.tolist())

            if g in self.uniprot_will_not_map:
                self.assertTrue(g in self.my_gempro.missing_uniprot_mapping)
                self.assertFalse(op.exists(op.join(gene_structure_folder, '{}.fasta'.format(u))))
                self.assertFalse(op.exists(op.join(gene_structure_folder, '{}.txt'.format(u))))
            else:
                uniprots = self.my_gempro.df_uniprot_metadata[self.my_gempro.df_uniprot_metadata.gene == g].uniprot.unique().tolist()
                for u in uniprots:
                    self.assertTrue(op.exists(op.join(gene_structure_folder, '{}.fasta'.format(u))))
                    self.assertTrue(op.exists(op.join(gene_structure_folder, '{}.txt'.format(u))))
                    # 3.
                    # Test that the.annotation['sequence'] attribute has UniProt information in it
                    u_prop = gene.protein.sequences.get_by_id(u)
                    self.assertTrue(u_prop.seq_len > 0)

    def test_x_map_uniprot_to_pdb(self):
        self.my_gempro.uniprot_mapping_and_metadata(model_gene_source='ENSEMBLGENOME_ID')
        self.my_gempro.set_representative_sequence()
        self.my_gempro.map_uniprot_to_pdb(seq_ident_cutoff=.9)

        look_at_this_gene = self.my_gempro.genes.get_by_id('b3734').protein
        self.assertTrue(look_at_this_gene.structures.has_id('3oaa'))
        look_at_this_structure = look_at_this_gene.structures.get_by_id('3oaa')
        # TODO: check parsing of taxonomy name
        # self.assertEqual(536056, look_at_this_structure.taxonomy_name)
        look_at_this_chain = look_at_this_structure.chains.get_by_id('A')
        self.assertEqual(513, look_at_this_chain.unp_end)
        self.assertEqual(1, look_at_this_chain.coverage)

    def test_x_blast_pdb(self):
        self.my_gempro.uniprot_mapping_and_metadata(model_gene_source='ENSEMBLGENOME_ID')
        self.my_gempro.set_representative_sequence()
        self.my_gempro.blast_seqs_to_pdb(seq_ident_cutoff=.8, all_genes=True)

        look_at_this_gene = self.my_gempro.genes.get_by_id('b2296').protein
        self.assertTrue(look_at_this_gene.structures.has_id('3slc'))
        look_at_this_structure = look_at_this_gene.structures.get_by_id('3slc')
        look_at_this_chain = look_at_this_structure.chains.get_by_id('A')
        self.assertEqual(0.9425, look_at_this_chain.hit_percent_ident)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.base_dir)