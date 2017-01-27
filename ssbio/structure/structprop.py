from ssbio.core.object import Object
from ssbio.structure.chainprop import ChainProp
import ssbio.sequence.utils.alignment
import os.path as op
from ssbio.structure.utils.pdbioext import PDBIOExt
from cobra.core import DictList
import ssbio.structure.properties
import ssbio.utils
import logging
import nglview as nv
import seaborn as sns
log = logging.getLogger(__name__)


class StructProp(Object):
    """Class for protein structural properties"""

    def __init__(self, ident, description=None, chains=None, mapped_chains=None, structure_file=None, file_type=None,
                 reference_seq=None, representative_chain=None, is_experimental=False, parse=False):
        Object.__init__(self, id=ident, description=description)

        self.is_experimental = is_experimental
        self.structure = None

        self.reference_seq = reference_seq
        self.reference_seq_top_coverage = 0

        # chains is a DictList of ChainProp objects
        self.chains = DictList()
        if chains:
            self.add_chain_ids(chains)

        # representative_chain is a pointer to the chain in self.chains that matches reference_seq
        self.representative_chain = None
        if representative_chain:
            self.representative_chain = self.chains.get_by_id(representative_chain)
        # mapped_chains is an ordered list of mapped chain IDs which would come from BLAST or the best_structures API
        self.mapped_chains = ssbio.utils.force_list(mapped_chains)

        self.file_type = file_type
        self.structure_path = structure_file
        if structure_file:
            self.load_structure_file(structure_file, file_type, parse)

    @property
    def structure_file(self):
        if not self.structure_path:
            return None
        return op.basename(self.structure_path)

    def load_structure_file(self, structure_file, file_type, parse=False):
        """Load a structure file and provide pointers to its location

        Args:
            structure_file: Path to structure file
            file_type: Type of structure file
            parse (bool): If the 3D coordinates should be parsed using Biopython
            get_chain_sequences (bool): If chain sequences should be parsed and stored in their corresponding ChainProp object

        """
        self.file_type = file_type
        self.structure_path = structure_file

        if parse:
            my_structure = self.parse_structure()
            if not my_structure:
                log.error('{}: unable to load structure file'.format(self.id))
                return
            log.debug('{}: parsed 3D coordinates into Biopython structure object'.format(self.id))

    def parse_structure(self):
        """Read the 3D coordinates of a structure file and save it in the structure attribute.

        Also create ChainProp objects in the chains attribute

        Args:
            get_chain_sequences (bool): If chain sequences should be parsed and stored in their corresponding ChainProp object

        Returns:
            Structure: Biopython structure object

        """
        if not self.structure_path:
            log.error('{}: no structure loaded, unable to parse'.format(self.id))
            return None
        else:
            # Add Biopython structure object
            self.structure = PDBIOExt(self.structure_path, file_type=self.file_type)

            # Add all chains to self.chains as ChainProp objects
            structure_chains = [x.id for x in self.structure.first_model.child_list]
            self.add_chain_ids(structure_chains)
            self._get_structure_seqs()

            # Also add all chains to self.mapped_chains ONLY if there are none specified
            if not self.mapped_chains:
                self.add_mapped_chain_ids(structure_chains)

        return self.structure

    def add_mapped_chain_ids(self, mapped_chains):
        """Add chains by ID into the mapped_chains attribute

        Args:
            mapped_chains (str, list): Chain ID or list of IDs

        """
        mapped_chains = ssbio.utils.force_list(mapped_chains)

        for c in mapped_chains:
            if c not in self.mapped_chains:
                self.mapped_chains.append(c)
                log.debug('{}: added to list of mapped chains'.format(c))
            else:
                log.debug('{}: chain already in list of mapped chains, not adding'.format(c))

    def add_chain_ids(self, chains):
        """Add chains by ID into the chains attribute

        Args:
            chains (str, list): Chain ID or list of IDs

        """
        chains = ssbio.utils.force_list(chains)

        for c in chains:
            # Check if chain ID is not in the structure
            if self.structure:
                if c not in self.structure.first_model:
                    raise ValueError('{}: specified chain not contained in model'.format(c))

            if self.chains.has_id(c):
                log.debug('{}: chain already present'.format(c))
            else:
                chain_prop = ChainProp(ident=c, pdb_parent=self.id)
                self.chains.append(chain_prop)
                log.debug('{}: added to chains list'.format(c))

    def _get_structure_seqs(self):
        """Store chain sequences in the corresponding ChainProp objects in the chains attribute

        Returns:
            DictList: All chain sequences as a DictList of SeqRecords

        """

        # Returns the structures sequences with Xs added
        structure_seqs = ssbio.structure.properties.residues.get_structure_seqrecords(self.structure.first_model)
        log.debug('{}: gathered chain sequences'.format(self.id))

        # Associate with ChainProps
        for seq_record in structure_seqs:
            log.debug('{}: adding chain sequence to ChainProp'.format(seq_record.id))
            my_chain = self.chains.get_by_id(seq_record.id)
            my_chain.seq_record = seq_record

    def align_reference_seq_to_mapped_chains(self, outdir=None, engine='needle', parse=False, force_rerun=False,
                                             **kwargs):
        """Run and store alignments of the reference sequence to chains in mapped_chains

        Alignments are stored in the reference_seq.structure_alignments attribute.

        """
        if not self.reference_seq:
            raise ValueError('{}: reference sequence not set'.format(self.id))

        # If no mapped chains have been specified and the structure hasn't been parsed, parse it and look at all chains
        if not self.mapped_chains and not self.structure:
            parsed = self.parse_structure()
            if not parsed:
                log.error('{}: unable to load structure file'.format(self.id))
                return

        for chain_id in self.mapped_chains:
            structure_id = '{}-{}'.format(self.id, chain_id)
            aln_id = '{}_{}'.format(self.reference_seq.id, structure_id)
            outfile = '{}.needle'.format(aln_id)

            if self.reference_seq.structure_alignments.has_id(aln_id):
                log.debug('{}: alignment already completed'.format(chain_id))
                continue

            log.debug('{}: aligning to reference sequence {}'.format(structure_id, self.reference_seq.id))

            chain_prop = self.chains.get_by_id(chain_id)
            if not chain_prop:
                raise ValueError('{}: chain sequence not parsed ')
            chain_seq_record = chain_prop.seq_record

            aln = ssbio.sequence.utils.alignment.pairwise_sequence_alignment(a_seq=self.reference_seq.seq_str,
                                                                             a_seq_id=self.reference_seq.id,
                                                                             b_seq=chain_seq_record,
                                                                             b_seq_id=structure_id,
                                                                             engine=engine,
                                                                             outdir=outdir,
                                                                             outfile=outfile,
                                                                             force_rerun=force_rerun)

            # Add an identifier to the MultipleSeqAlignment object for storage in a DictList
            aln.id = aln_id
            aln.annotations['a_seq'] = self.reference_seq.id
            aln.annotations['b_seq'] = structure_id
            aln.annotations['structure_id'] = self.id
            aln.annotations['chain_id'] = chain_id

            if parse:
                aln_df = ssbio.sequence.utils.alignment.get_alignment_df(a_aln_seq=str(list(aln)[0].seq),
                                                                         b_aln_seq=str(list(aln)[1].seq))
                aln.annotations['mutations'] = ssbio.sequence.utils.alignment.get_mutations(aln_df)
                aln.annotations['deletions'] = ssbio.sequence.utils.alignment.get_deletions(aln_df)
                aln.annotations['insertions'] = ssbio.sequence.utils.alignment.get_insertions(aln_df)

            self.reference_seq.structure_alignments.append(aln)

    def sequence_quality_checker(self, seq_ident_cutoff=0.5, allow_missing_on_termini=0.2,
                                 allow_mutants=True, allow_deletions=False,
                                 allow_insertions=False, allow_unresolved=True):
        """Set the representative chain based on sequence quality checks to the reference sequence.

        Args:
            seq_ident_cutoff:
            allow_missing_on_termini:
            allow_mutants:
            allow_deletions:
            allow_insertions:
            allow_unresolved:

        Returns:

        """
        for alignment in self.reference_seq.structure_alignments:
            chain_id = alignment.annotations['chain_id']

            # Compare representative sequence to structure sequence using the alignment
            found_good_chain = ssbio.structure.properties.quality.sequence_checker(reference_seq_aln=alignment[0],
                                                                        structure_seq_aln=alignment[1],
                                                                        seq_ident_cutoff=seq_ident_cutoff,
                                                                        allow_missing_on_termini=allow_missing_on_termini,
                                                                        allow_mutants=allow_mutants,
                                                                        allow_deletions=allow_deletions,
                                                                        allow_insertions=allow_insertions,
                                                                        allow_unresolved=allow_unresolved)

            # If found_good_pdb = True, set as representative chain
            # If not, move on to the next potential chain
            if found_good_chain:
                self.representative_chain = self.chains.get_by_id(chain_id)
                self.reference_seq_top_coverage = alignment.annotations['percent_identity']
                log.debug('{}: chain {} set as representative'.format(self.id, chain_id))
                return self.representative_chain
        else:
            log.debug('{}: no chains meet quality checks'.format(self.id))
            return None

    def get_dict_with_chain(self, chain, only_keys=None, chain_keys=None, exclude_attributes=None, df_format=False):
        """get_dict method which incorporates attributes found in a specific chain. Does not overwrite any attributes
            in the original StructProp.

        Args:
            chain:
            only_keys:
            chain_keys:
            exclude_attributes:
            df_format:

        Returns:
            dict: attributes of StructProp + the chain specified

        """

        # Choose attributes to return, return everything in the object if a list is not specified
        if not only_keys:
            keys = list(self.__dict__.keys())
        else:
            keys = ssbio.utils.force_list(only_keys)

        # Remove keys you don't want returned
        if exclude_attributes:
            exclude_attributes = ssbio.utils.force_list(exclude_attributes)
            for x in exclude_attributes:
                if x in keys:
                    keys.remove(x)
        else:
            exclude_attributes = []

        exclude_attributes.extend(['mapped_chains', 'chains'])

        final_dict = {k: v for k, v in Object.get_dict(self, only_keys=keys, exclude_attributes=exclude_attributes,
                                                       df_format=df_format).items()}


        chain_prop = self.chains.get_by_id(chain)
        # Filter out keys that show up in StructProp
        if not chain_keys:
            chain_keys = [x for x in chain_prop.get_dict().keys() if x not in final_dict]

        chain_dict = chain_prop.get_dict(only_keys=chain_keys, df_format=df_format)
        final_dict.update(chain_dict)

        return final_dict

    def get_disulfide_bridges(self, threshold=3.0):
        """Run Biopython's search_ss_bonds to find potential disulfide bridges for each chain and store in ChainProp.
        """
        if not self.structure:
            parsed = self.parse_structure()
            if not parsed:
                log.error('{}: unable to open structure to find S-S bridges'.format(self.id))
                return

        disulfide_bridges = ssbio.structure.properties.residues.search_ss_bonds(self.structure.first_model,
                                                                                threshold=threshold)
        if not disulfide_bridges:
            log.debug('{}: no disulfide bridges'.format(self.id))

        for chain, bridges in disulfide_bridges.items():
            self.representative_chain.seq_record.annotations['SSBOND-biopython'] = disulfide_bridges[self.representative_chain.id]
            log.debug('{}: found {} disulfide bridges'.format(chain, len(bridges)))
            log.debug('{}: stored disulfide bridges in seq_record letter_annotations'.format(chain))

    def get_residue_depths(self, outdir, force_rerun=False):
        """Run MSMS on this structure and store the residue depths/ca depths in the corresponding ChainProp SeqRecords
        """
        if not self.structure:
            parsed = self.parse_structure()
            if not parsed:
                log.error('{}: unable to open structure to run MSMS'.format(self.id))
                return

        log.debug('{}: running MSMS'.format(self.id))
        msms_results = ssbio.structure.properties.msms.get_msms_df(model=self.structure.first_model,
                                                        pdb_file=self.structure_path,
                                                        outdir=outdir, force_rerun=force_rerun)
        if msms_results.empty:
            log.error('{}: unable to run MSMS'.format(self.id))
            return

        chains = msms_results.chain.unique()

        for chain in chains:
            res_depths = msms_results[msms_results.chain == chain].res_depth.tolist()
            ca_depths = msms_results[msms_results.chain == chain].ca_depth.tolist()

            chain_prop = self.chains.get_by_id(chain)
            chain_seq = chain_prop.seq_record

            # Making sure the X's are filled in
            res_depths = ssbio.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                           new_seq=res_depths,
                                                                           fill_with=float('Inf'))

            ca_depths = ssbio.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                          new_seq=ca_depths,
                                                                          fill_with=float('Inf'))

            chain_prop.seq_record.letter_annotations['RES_DEPTH-msms'] = res_depths
            chain_prop.seq_record.letter_annotations['CA_DEPTH-msms'] = ca_depths
            log.debug('{}: stored residue depths in chain seq_record letter_annotations'.format(chain))

    def get_dssp_annotations(self, outdir, force_rerun=False):
        """Run DSSP on this structure and store the DSSP annotations in the corresponding ChainProp SeqRecords

        Args:
            outdir (str): Path to where DSSP dataframe will be stored.
            force_rerun (bool): If DSSP results should be recalculated

        """
        if not self.structure:
            parsed = self.parse_structure()
            if not parsed:
                log.error('{}: unable to open structure to run DSSP'.format(self.id))
                return

        log.debug('{}: running DSSP'.format(self.id))
        dssp_results = ssbio.structure.properties.dssp.get_dssp_df(model=self.structure.first_model,
                                                                   pdb_file=self.structure_path,
                                                                   outdir=outdir,
                                                                   force_rerun=force_rerun)

        if dssp_results.empty:
            log.error('{}: unable to run DSSP'.format(self.id))
            return

        chains = dssp_results.chain.unique()
        dssp_summary = ssbio.structure.properties.dssp.secondary_structure_summary(dssp_results)

        for chain in chains:
            ss = dssp_results[dssp_results.chain == chain].ss.tolist()
            exposure_rsa = dssp_results[dssp_results.chain == chain].exposure_rsa.tolist()
            exposure_asa = dssp_results[dssp_results.chain == chain].exposure_asa.tolist()
            phi = dssp_results[dssp_results.chain == chain].phi.tolist()
            psi = dssp_results[dssp_results.chain == chain].psi.tolist()

            chain_prop = self.chains.get_by_id(chain)
            chain_seq = chain_prop.seq_record

            # Making sure the X's are filled in
            ss = ssbio.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                              new_seq=ss,
                                                                              fill_with='-')
            exposure_rsa = ssbio.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                             new_seq=exposure_rsa,
                                                                             fill_with=float('Inf'))
            exposure_asa = ssbio.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                             new_seq=exposure_asa,
                                                                             fill_with=float('Inf'))
            phi = ssbio.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                    new_seq=phi,
                                                                    fill_with=float('Inf'))
            psi = ssbio.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                    new_seq=psi,
                                                                    fill_with=float('Inf'))

            chain_prop.seq_record.annotations.update(dssp_summary[chain])

            chain_prop.seq_record.letter_annotations['SS-dssp'] = ss
            chain_prop.seq_record.letter_annotations['RSA-dssp'] = exposure_rsa
            chain_prop.seq_record.letter_annotations['ASA-dssp'] = exposure_asa
            chain_prop.seq_record.letter_annotations['PHI-dssp'] = phi
            chain_prop.seq_record.letter_annotations['PSI-dssp'] = psi
            log.debug('{}: stored DSSP annotations in chain seq_record letter_annotations'.format(chain))

    def view_structure(self, opacity=1.0, gui=False):
        if not self.structure_path:
            raise ValueError("Structure file not loaded")
        view = nv.show_structure_file(self.structure_path, gui=gui)
        view.clear_representations()
        view.add_cartoon(selection='protein', color='silver', opacity=opacity)
        return view

    def view_structure_with_mutations(self, mutations, color='red', unique_colors=False, structure_opacity=0.5,
                                      opacity_range=(0.5,1), scale_range=(.7, 10), gui=False):
        """Input a list of residue numbers to view on the structure. Or input a dictionary of residue numbers to counts.

        Args:
            mutations (int, list)
            color:
            opacity_range:
            scale_range:

        Returns:

        """
        opacity_dict = ssbio.utils.scale_calculator(opacity_range[0], mutations, rescale=opacity_range)
        scale_dict = ssbio.utils.scale_calculator(scale_range[0], mutations, rescale=scale_range)

        view = self.view_structure(opacity=structure_opacity, gui=gui)

        if isinstance(mutations, list):
            unique_mutations = list(set(mutations))
        elif isinstance(mutations, dict):
            unique_mutations = list(mutations.keys())

        # TODO: add color by letter_annotations!

        colors = sns.color_palette("hls", len(unique_mutations)).as_hex()

        for i, x in enumerate(unique_mutations):
            if isinstance(x, tuple):
                to_show = ''
                for mut in x:
                    to_show += '{} or '.format(mut)
                to_show = to_show.strip(' or ')
            else:
                to_show = x

            if unique_colors:
                view.add_ball_and_stick(selection='not hydrogen and {}'.format(to_show),
                                        color=colors[i], opacity=opacity_dict[x], scale=scale_dict[x])
            else:
                view.add_ball_and_stick(selection='not hydrogen and {}'.format(to_show),
                                        color=color, opacity=opacity_dict[x], scale=scale_dict[x])

        return view

