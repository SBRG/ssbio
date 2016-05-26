import os.path as op

class GEMProPrep():
    '''
    Prepares folders for GEM-PRO creation
    '''

    def __init__(self, gem_name, root_dir):
        self.root_dir = root_dir

        self.model_dir = op.join(root_dir, gem_name)

        # data_frames - directory where all data frames will be stored (all stages)
        self.data_frames = op.join(model_dir, 'data_frames')

        # figures - directory where all figures will be stored (all stages)
        figures = op.join(model_dir, 'figures')

        # model_files - directory where original gems and gem-related files are stored
        model_files = op.join(model_dir, 'model_files')

        # structure_files - directory where structure related files will be downloaded/are located
        struct_files = op.join(model_dir, 'structure_files')
        struct_exp_files = op.join(struct_files, 'experimental')
        struct_homology_files = op.join(struct_files, 'homology_models')

        seq_files = op.join(model_dir, 'sequence_files')
        seq_homology_files = op.join(seq_files, 'homology_model_sequences')
        seq_uniprot_files = op.join(seq_files, 'uniprot_sequences')
        seq_uniprot_metadata = op.join(seq_files, 'uniprot_metadata')
        seq_align_files = op.join(seq_files, 'alignment')
        seq_pdb_files = op.join(seq_files, 'pdb_sequences')

    def prep_folders(self, seq):
        """Take in a sequence string and prepares the folder for it to run ITASSER
        """
        itasser_dir = op.join(self.root_dir, self.ident)