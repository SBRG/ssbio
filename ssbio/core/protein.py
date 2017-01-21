import os
import pandas as pd
from cobra.core import DictList
import ssbio.databases.kegg
import ssbio.databases.uniprot
from ssbio.sequence import SeqProp
from ssbio.databases.kegg import KEGGProp
from ssbio.databases.uniprot import UniProtProp
from ssbio.structure import StructProp
from ssbio.databases.pdb import PDBProp
from ssbio.structure.homology.itasser.itasserprop import ITASSERProp
import ssbio.sequence.utils.fasta
import ssbio.sequence.utils.alignment
import ssbio.utils
import requests
from ssbio.structure.utils.cleanpdb import CleanPDB
import ssbio.structure.properties.quality
import numpy as np
import logging
log = logging.getLogger(__name__)


class Protein(object):
    """Basic definition of a protein"""

    _representative_sequence_attributes = ['gene', 'uniprot', 'kegg', 'pdbs', 'sequence_len',
                                           'sequence_file', 'sequence_path', 'metadata_file', 'metadata_path']
    _representative_structure_attributes = ['is_experimental', 'reference_seq_top_coverage', 'date', 'description',
                                            'resolution','taxonomy_name']

    def __init__(self, ident):
        self.id = ident
        self.sequences = DictList()
        self.structures = DictList()
        self.representative_sequence = None
        self.representative_structure = None

        # TODO: define this instead of doing op.join(gempro.sequence_dir, gene.id) for a lot of things
        # self.sequence_dir = seq_dir
        # self.structure_dir = struct_dir

    @property
    def num_structures(self):
        return len(self.structures)

    @property
    def num_structures_experimental(self):
        return len(self.get_experimental_structures())

    @property
    def num_structures_homology(self):
        return len(self.get_homology_models())

    def get_experimental_structures(self):
        # Return a DictList of all experimental structures in self.structures
        return DictList(x for x in self.structures if x.is_experimental)

    def get_homology_models(self):
        # Return a DictList of all homology models in self.structures
        return DictList(x for x in self.structures if not x.is_experimental)

    def filter_sequences(self, seq_type):
        """Get a DictList of only specified mappings in the sequences attribute.

        Args:
            seq_type: Object type

        Returns:
            DictList: of Object type mappings only

        """
        return DictList(x for x in self.sequences if isinstance(x, seq_type))

    def load_kegg(self, kegg_id, kegg_organism_code=None, kegg_seq_file=None, kegg_metadata_file=None,
                  set_as_representative=False,
                  download=False, outdir=None, force_rerun=False):
        """Load existing KEGG ID, sequence, and metadata files into the sequences attribute.

        Args:
            kegg_id:
            kegg_organism_code: KEGG organism code to prepend to the kegg_id if not part of it already.
                Example: "eco:b1244", eco is the organism code
            kegg_seq_file:
            kegg_metadata_file:
            set_as_representative:

        Returns:
            KEGGProp: object contained in the sequences attribute

        """
        if kegg_organism_code:
            kegg_id = kegg_organism_code + ':' + kegg_id

        # If KEGG ID has already been mapped, just return the KEGGProp object in sequences
        if kegg_id not in self.sequences:
            kegg_prop = KEGGProp(kegg_id, kegg_seq_file, kegg_metadata_file)
            if download:
                kegg_prop.download_seq_file(outdir, force_rerun)
                kegg_prop.download_metadata_file(outdir, force_rerun)

            if kegg_prop.sequence_path:
                # Check if KEGG sequence matches a potentially set representative sequence
                # Do not add any info if a UniProt ID was already mapped though, we want to use that
                if self.representative_sequence:
                    if not self.representative_sequence.uniprot:
                        if kegg_prop.equal_to(self.representative_sequence):
                            # Update the representative sequence field with KEGG metadata
                            self.representative_sequence.update(kegg_prop.get_dict())
                        else:
                            # TODO: add option to use manual or kegg sequence if things do not match
                            log.warning('{}: representative sequence does not match mapped KEGG sequence.'.format(self.id))

            if set_as_representative:
                self._representative_sequence_setter(kegg_prop)

            self.sequences.append(kegg_prop)

        return self.sequences.get_by_id(kegg_id)

    def load_uniprot(self, uniprot_id, uniprot_seq_file=None, uniprot_metadata_file=None,
                     set_as_representative=False,
                     download=False, outdir=None, force_rerun=False):
        """Load a UniProt ID and associated sequence/metadata files into the sequences attribute.

        Args:
            uniprot_id:
            uniprot_seq_file:
            uniprot_metadata_file:

        """
        if uniprot_id not in self.sequences:
            uniprot_prop = UniProtProp(uniprot_id, uniprot_seq_file, uniprot_metadata_file)
            if download:
                uniprot_prop.download_seq_file(outdir, force_rerun)
                uniprot_prop.download_metadata_file(outdir, force_rerun)

            # Also check if UniProt sequence matches a potentially set representative sequence
            if self.representative_sequence:
                # Test equality
                if uniprot_prop.equal_to(self.representative_sequence):
                    # Update the representative sequence field with UniProt metadata
                    self.representative_sequence.update(uniprot_prop.get_dict())
                else:
                    # TODO: add option to use manual or uniprot sequence if things do not match
                    log.warning('{}: representative sequence does not match mapped UniProt sequence'.format(self.id))

            if set_as_representative:
                self._representative_sequence_setter(uniprot_prop)

            self.sequences.append(uniprot_prop)
        return self.sequences.get_by_id(uniprot_id)

    def load_manual_sequence_file(self, ident, seq_file, set_as_representative=False):
        """Load a manual sequence given as a FASTA file and optionally set it as the representative sequence.
            Also store it in the sequences attribute.

        Args:
            ident:
            seq_file:
            set_as_representative:

        """
        manual_sequence = SeqProp(ident=ident, sequence_file=seq_file)
        self.sequences.append(manual_sequence)

        if set_as_representative:
            self._representative_sequence_setter(manual_sequence)

        return self.sequences.get_by_id(ident)

    def load_manual_sequence_str(self, ident, seq_str, outdir=None, set_as_representative=False):
        """Load a manual sequence given as a string and optionally set it as the representative sequence.

        Also store it in the sequences attribute.

        Args:
            ident:
            seq_str:
            set_as_representative:

        """
        seq_file = ssbio.sequence.utils.fasta.write_fasta_file(indict={ident: seq_str}, outname=ident, outdir=outdir)
        manual_sequence = SeqProp(ident=ident, sequence_file=seq_file)
        self.sequences.append(manual_sequence)

        if set_as_representative:
            self._representative_sequence_setter(manual_sequence)

        return self.sequences.get_by_id(ident)

    def _representative_sequence_setter(self, seq_prop):
        """Make a copy of a SeqProp object and store it as the representative. Only keep certain attributes"""
        self.representative_sequence = SeqProp(ident=seq_prop.id)
        self.representative_sequence.update(seq_prop.get_dict(), only_keys=self._representative_sequence_attributes)

    def set_representative_sequence(self):
        """Consolidate sequences and set a single representative sequence.
        """

        if len(self.sequences) == 0:
            log.error('{}: no sequences mapped'.format(self.id))
            return self.representative_sequence

        kegg_mappings = self.filter_sequences(KEGGProp)
        if len(kegg_mappings) > 0:
            kegg_to_use = kegg_mappings[0]
            if len(kegg_mappings) > 1:
                log.warning('{}: multiple KEGG mappings found, using the first entry {}'.format(self.id, kegg_to_use.id))

        uniprot_mappings = self.filter_sequences(UniProtProp)

        # If a representative sequence has already been set, nothing needs to be done
        if self.representative_sequence:
            log.debug('{}: representative sequence already set'.format(self.id))

        # If there is a KEGG annotation and no UniProt annotations, set KEGG as representative
        elif len(kegg_mappings) > 0 and len(uniprot_mappings) == 0:
            # TODO: double check if everything is set properly in seqprop
            self._representative_sequence_setter(kegg_to_use)
            log.debug('{}: representative sequence set from KEGG ID {}'.format(self.id, kegg_to_use.id))

        # If there are UniProt annotations and no KEGG annotations, set UniProt as representative
        elif len(kegg_mappings) == 0 and len(uniprot_mappings) > 0:
            # If there are multiple uniprots rank them by the sum of reviewed (bool) + num_pdbs
            # This way, UniProts with PDBs get ranked to the top, or if no PDBs, reviewed entries
            u_ranker = []
            for u in uniprot_mappings:
                u_ranker.append((u.id, u.ranking_score()))
            sorted_by_second = sorted(u_ranker, key=lambda tup: tup[1], reverse=True)
            best_u_id = sorted_by_second[0][0]

            best_u = uniprot_mappings.get_by_id(best_u_id)
            self._representative_sequence_setter(best_u)
            log.debug('{}: Representative sequence set from UniProt ID {}'.format(self.id, best_u_id))

        # If there are both UniProt and KEGG annotations...
        elif len(kegg_mappings) > 0 and len(uniprot_mappings) > 0:
            # Use KEGG if the mapped UniProt is unique, and it has PDBs
            if kegg_to_use.num_pdbs() > 0 and not uniprot_mappings.has_id(kegg_to_use.uniprot):
                self._representative_sequence_setter(kegg_to_use)
                log.debug('{}: Representative sequence set from KEGG ID {}'.format(self.id, kegg_to_use.id))
            else:
                # If there are multiple uniprots rank them by the sum of reviewed (bool) + num_pdbs
                u_ranker = []
                for u in uniprot_mappings:
                    u_ranker.append((u.id, u.ranking_score()))
                sorted_by_second = sorted(u_ranker, key=lambda tup: tup[1], reverse=True)
                best_u_id = sorted_by_second[0][0]

                best_u = uniprot_mappings.get_by_id(best_u_id)
                self._representative_sequence_setter(best_u)
                log.debug('{}: Representative sequence set from UniProt ID {}'.format(self.id, best_u_id))

        return self.representative_sequence

    def align_sequences_to_representative(self, outdir=None, engine='needle', parse=False, force_rerun=False, **kwargs):
        """Align all sequences in the sequences attribute to the representative sequence.

        Stores the alignments the representative_sequence.sequence_alignments DictList

        Args:
            outfile:
            outdir:
            engine:
            parse: Store locations of mutations, insertions, and deletions in the alignment object (as an annotation)
            force_rerun:
            **kwargs: Other options for sequence alignment (gap penalties, etc) See:
                ssbio.sequence.utils.alignment.pairwise_sequence_alignment

        """
        for seq in self.sequences:
            aln_id = '{}_{}'.format(self.id, seq.id)
            outfile = '{}.needle'.format(aln_id)

            if self.representative_sequence.sequence_alignments.has_id(aln_id):
                log.debug('{}: alignment already completed'.format(seq.id))
                continue

            aln = ssbio.sequence.utils.alignment.pairwise_sequence_alignment(a_seq=self.representative_sequence.get_seq_str(),
                                                                             a_seq_id=self.id,
                                                                             b_seq=seq.get_seq_str(),
                                                                             b_seq_id=seq.id,
                                                                             engine=engine,
                                                                             outdir=outdir,
                                                                             outfile=outfile,
                                                                             force_rerun=force_rerun)
            # Add an identifier to the MultipleSeqAlignment object for storage in a DictList
            aln.id = aln_id
            aln.annotations['a_seq'] = self.representative_sequence.id
            aln.annotations['b_seq'] = seq.id

            if parse:
                aln_df = ssbio.sequence.utils.alignment.get_alignment_df(a_aln_seq=str(list(aln)[0].seq),
                                                                         b_aln_seq=str(list(aln)[1].seq))
                aln.annotations['mutations'] = ssbio.sequence.utils.alignment.get_mutations(aln_df)
                aln.annotations['deletions'] = ssbio.sequence.utils.alignment.get_deletions(aln_df)
                aln.annotations['insertions'] = ssbio.sequence.utils.alignment.get_insertions(aln_df)

            self.representative_sequence.sequence_alignments.append(aln)

    def load_pdb(self, pdb_id, mapped_chains=None, pdb_file=None, file_type=None, set_as_representative=False, parse=False):
        """Load a PDB ID into the structures attribute.

        Args:
            pdb_id (str): PDB ID
            mapped_chains (str, list): Chain ID or list of IDs which you are interested in
            pdb_file (str): Path to PDB file
            file_type (str): Type of PDB file
            set_as_representative (bool): If this structure should be set as the representative structure
            parse (bool): If the structure's 3D coordinates and chains should be parsed

        Returns:
            PDBProp: The object that is now contained in the structures attribute

        """
        pdb_id = pdb_id.lower()

        if self.structures.has_id(pdb_id):
            log.debug('{}: PDB ID already present in list of structures'.format(pdb_id))
            pdb = self.structures.get_by_id(pdb_id)

            if mapped_chains:
                pdb.add_mapped_chain_ids(mapped_chains)
            if pdb_file:
                pdb.load_structure_file(pdb_file, file_type, parse)
        else:
            pdb = PDBProp(ident=pdb_id, chains=mapped_chains, structure_file=pdb_file, file_type=file_type, parse=parse,
                          reference_seq=self.representative_sequence)
            if mapped_chains:
                pdb.add_mapped_chain_ids(mapped_chains)
            self.structures.append(pdb)

        if set_as_representative:
            self._representative_structure_setter(pdb)

        return self.structures.get_by_id(pdb_id)

    def load_itasser_folder(self, ident, itasser_folder, set_as_representative=False, create_dfs=False, parse=False):
        itasser = ITASSERProp(ident, itasser_folder, create_dfs=create_dfs, parse=parse,
                              reference_seq=self.representative_sequence)
        if self.structures.has_id(itasser.id):
            log.warning('{}: already present in list of structures'.format(itasser.id))
            self.structures.get_by_id(ident)
        else:
            self.structures.append(itasser)

        if set_as_representative:
            self._representative_structure_setter(itasser)

        return self.structures.get_by_id(ident)

    def load_homology_model(self, ident, pdb_file=None, set_as_representative=False, parse=False):
        pdb = StructProp(ident=ident, structure_file=pdb_file, is_experimental=False, parse=parse,
                         reference_seq=self.representative_sequence)
        self.structures.append(pdb)

        if set_as_representative:
            self._representative_structure_setter(pdb)

        return self.structures.get_by_id(ident)

    def _representative_structure_setter(self, struct_prop, keep_chain, new_id=None, clean=True, out_suffix='_clean', outdir=None):
        """Set the representative structure by cleaning it and copying over attributes of the original structure.

        Args:
            struct_prop (StructProp): StructProp object to set as representative
            keep_chain (str, list): List of chains to keep
            new_id (str): New ID to call this structure, for example 1abc-D to represent PDB 1abc, chain D
            clean (bool): If the PDB file should be cleaned (see ssbio.structure.utils.cleanpdb)
            out_suffix (str): Suffix to append to clean PDB file
            outdir (str): Path to output directory

        Returns:
            StructProp: representative structure

        """
        # TODO: which attributes to keep?
        if not outdir:
            outdir = os.getcwd()

        if not new_id:
            new_id = struct_prop.id

        # Parse the structure if it hasn't already been, so chains will be added
        if not struct_prop.structure:
            parsed = struct_prop.parse_structure()
            if not parsed:
                log.error('{}: can\'t parse, unable to set as representative structure'.format(self.id))
                return

        # If the structure is to be cleaned, and which chain to keep
        if clean:
            custom_clean = CleanPDB(keep_chains=keep_chain)

            final_pdb = struct_prop.structure.write_pdb(custom_selection=custom_clean,
                                                        out_suffix=out_suffix,
                                                        out_dir=outdir)
        else:
            final_pdb = struct_prop.structure_path

        self.representative_structure = StructProp(ident=new_id, structure_file=final_pdb,
                                                   reference_seq=self.representative_sequence,
                                                   chains=keep_chain, representative_chain=keep_chain, mapped_chains=keep_chain,
                                                   file_type='pdb', parse=True)
        # structure needs to be excluded when using get_dict because of too much recursion
        self.representative_structure.update(struct_prop.get_dict_with_chain(chain=keep_chain,
                                                                             exclude_attributes='structure'),
                                             only_keys=self._representative_structure_attributes,
                                             overwrite=True)

        # STORE REPRESENTATIVE CHAIN RESNUMS in the representative sequence seqrecord letter_annotations
        # Get the alignment
        alnid = '{}_{}'.format(self.representative_sequence.id, self.representative_structure.id)
        aln = self.representative_structure.reference_seq.structure_alignments.get_by_id(alnid)
        # Get the mapping and store it in .seq_record.letter_annotations['repchain_resnums']
        aln_df = ssbio.sequence.utils.alignment.get_alignment_df(aln[0], aln[1])
        repchain_resnums = aln_df[pd.notnull(aln_df.id_a_pos)].id_b_pos.tolist()
        self.representative_sequence.seq_record.letter_annotations['repchain_resnums'] = repchain_resnums

    def set_representative_structure(self, seq_outdir, struct_outdir, engine='needle', seq_ident_cutoff=0.5,
                                     always_use_homology=False,
                                     allow_missing_on_termini=0.2, allow_mutants=True, allow_deletions=False,
                                     allow_insertions=False, allow_unresolved=True, force_rerun=False):
        """Set a representative structure from the structures in self.structures

        Args:
            always_use_homology:
            sort_homology_by:
            allow_missing_on_termini:
            allow_mutants:
            allow_deletions:
            allow_insertions:
            allow_unresolved:
            force_rerun:

        Returns:

        """
        log.debug('{}: setting representative structure'.format(self.id))

        if len(self.structures) == 0:
            log.debug('{}: no structures available'.format(self.id))
            return None

        if not self.representative_sequence:
            log.error('{}: no representative sequence to compare structures to'.format(self.id))
            return None

        if self.representative_structure and not force_rerun:
            log.debug('{}: representative structure already set'.format(self.id))
            return self.representative_structure

        has_homology = False
        has_pdb = False
        use_homology = False
        use_pdb = False

        if self.num_structures_homology > 0:
            has_homology = True
        if self.num_structures_experimental > 0:
            has_pdb = True

        # If we mark to always use homology, use it if it exists
        if always_use_homology:
            if has_homology:
                use_homology = True
            elif has_pdb:
                use_pdb = True
        # If we don't always want to use homology, use PDB if it exists
        else:
            if has_homology and has_pdb:
                use_pdb = True
                use_homology = True
            elif has_homology and not has_pdb:
                use_homology = True
            elif has_pdb and not has_homology:
                use_pdb = True

        if use_pdb:
            # Put PDBs through QC/QA
            log.debug('{}: checking quality of experimental structures'.format(self.id))
            all_pdbs = self.get_experimental_structures()

            for pdb in all_pdbs:
                # Download the structure and parse it
                # This will add all chains to the mapped_chains attribute if there are none
                try:
                    pdb.download_structure_file(outdir=struct_outdir, force_rerun=force_rerun, parse=True)
                    # TODO: add global flag of PDB file type and adjust for downloading the header here (and other places where pdb is downloaded)
                    # pdb.download_cif_header_file(outdir=struct_outdir)
                except requests.exceptions.HTTPError:
                    log.error('{}: structure file could not be downloaded'.format(pdb))
                    continue
                pdb.align_reference_seq_to_mapped_chains(outdir=seq_outdir, engine=engine, parse=False, force_rerun=force_rerun)
                best_chain = pdb.sequence_quality_checker(seq_ident_cutoff=seq_ident_cutoff,
                                                          allow_missing_on_termini=allow_missing_on_termini,
                                                          allow_mutants=allow_mutants, allow_deletions=allow_deletions,
                                                          allow_insertions=allow_insertions, allow_unresolved=allow_unresolved)

                if best_chain:
                    self._representative_structure_setter(struct_prop=pdb,
                                                          new_id='{}-{}'.format(pdb.id, best_chain.id),
                                                          clean=True,
                                                          out_suffix='-{}_clean'.format(best_chain.id),
                                                          keep_chain=best_chain.id,
                                                          outdir=struct_outdir)
                    log.debug('{}-{}: set as representative structure'.format(pdb.id, best_chain.id))
                    return self.representative_structure
            else:
                log.debug('{}: no experimental structures meet cutoffs'.format(self.id))

        # If we are to use homology, save its information in the representative structure field
        if use_homology:
            log.debug('{}: checking quality of homology models'.format(self.id))
            all_models = self.get_homology_models()

            # TODO: homology models are not ordered in any other way other than how they are loaded,
            # rethink this for multiple homology models
            for homology in all_models:
                if not homology.structure_path:
                    log.debug('{}: no homology structure file'.format(self.id))
                    continue

                homology.align_reference_seq_to_mapped_chains(outdir=seq_outdir, engine=engine, parse=False,
                                                              force_rerun=force_rerun)
                best_chain = homology.sequence_quality_checker(seq_ident_cutoff=seq_ident_cutoff,
                                                               allow_missing_on_termini=allow_missing_on_termini,
                                                               allow_mutants=allow_mutants,
                                                               allow_deletions=allow_deletions,
                                                               allow_insertions=allow_insertions,
                                                               allow_unresolved=allow_unresolved)

                if best_chain:
                    self._representative_structure_setter(struct_prop=homology,
                                                          new_id='{}-{}'.format(homology.id, best_chain.id),
                                                          clean=True,
                                                          out_suffix='-{}_clean'.format(best_chain.id),
                                                          keep_chain=best_chain.id,
                                                          outdir=struct_outdir)
                    log.debug('{}-{}: set as representative structure'.format(homology.id, best_chain.id))
                    return self.representative_structure

        else:
            log.warning('{}: no representative structure found'.format(self.id))
            return None

    def map_repseq_resnums_to_repchain_index(self, resnums):
        """Map a residue number in the representative_sequence to an index in the representative_chain

        Use this to get the indices of the repchain to get structural properties at a specific residue number.

        Args:
            resnums (int, list): Residue numbers in the representative sequence

        Returns:
            dict: Mapping of resnums to indices

        """
        resnums = ssbio.utils.force_list(resnums)

        repchain_resnum_mapping = self.representative_sequence.seq_record.letter_annotations['repchain_resnums']

        to_repchain_index = {}
        for x in resnums:
            ix = repchain_resnum_mapping[x - 1] - 1

            if np.isnan(ix):
                log.debug('{}, {}: no equivalent residue found in structure'.format(self.id, x))
            else:
                to_repchain_index[x] = int(ix)

        return to_repchain_index

    def map_repseq_resnums_to_structure_resnums(self, resnums):
        """Map a residue number in the representative_sequence to the actual structure file's residue number

        Args:
            resnums (int, list): Residue numbers in the representative sequence

        Returns:
            dict: Mapping of resnums to structure residue IDs

        """
        resnums = ssbio.utils.force_list(resnums)

        mapping_to_repchain_index = self.map_repseq_resnums_to_repchain_index(resnums)
        repchain_structure_mapping = self.representative_structure.representative_chain.seq_record.letter_annotations['structure_resnums']

        to_structure_resnums = {}
        for k, v in mapping_to_repchain_index.items():
            rn = repchain_structure_mapping[v]

            if rn[1] == float('Inf'):
                log.debug('{}, {}: structure file does not contain coordinates for this residue'.format(self.id, k))
            else:
                to_structure_resnums[k] = rn

        return to_structure_resnums

    def view_all_mutations(self, grouped=False, color='red', unique_colors=True, structure_opacity=0.5,
                           opacity_range=(0.8,1), scale_range=(1,5), gui=False):
        """Map all sequence alignment mutations to the structure.

        Args:
            grouped (bool): If groups of mutations should be colored and sized together
            color (str): Color of the mutations (overridden if unique_colors=True)
            unique_colors (bool): If each mutation/mutation group should be colored uniquely
            opacity_range (tuple): Min/max opacity values (mutations that show up more will be opaque)
            scale_range (tuple): Min/max size values (mutations that show up more will be bigger)

        Returns:
            NGLviewer object

        """
        single, fingerprint = self.representative_sequence.sequence_mutation_summary()

        single_lens = {k: len(v) for k, v in single.items()}
        single_map_to_structure = {}
        for k, v in single_lens.items():
            resnum = int(k[1])
            resnum_to_structure = self.map_repseq_resnums_to_structure_resnums(resnum)
            new_key = resnum_to_structure[resnum][1]
            single_map_to_structure[new_key] = v

        if not grouped:
            view = self.representative_structure.view_structure_with_mutations(single_map_to_structure,
                                                                               color=color, unique_colors=unique_colors,
                                                                               structure_opacity=structure_opacity,
                                                                               opacity_range=opacity_range,
                                                                               scale_range=scale_range,
                                                                               gui=gui)
            return view

        else:
            fingerprint_lens = {k: len(v) for k, v in fingerprint.items()}
            fingerprint_map_to_structure = {}
            for k, v in fingerprint_lens.items():
                k_list = [int(x[1]) for x in k]
                resnums_to_structure = self.map_repseq_resnums_to_structure_resnums(k_list)
                new_key = tuple(y[1] for y in resnums_to_structure.values())
                fingerprint_map_to_structure[new_key] = v

            view = self.representative_structure.view_structure_with_mutations(fingerprint_map_to_structure,
                                                                               color=color, unique_colors=unique_colors,
                                                                               opacity_range=opacity_range,
                                                                               scale_range=scale_range,
                                                                               gui=gui)
            return view
