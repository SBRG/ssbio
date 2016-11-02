import glob
import os

import numpy as np
import pandas as pd
import ssbio.sequence.alignment
from ssbio.utils import percentage_to_float

import logging
log = logging.getLogger(__name__)

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

def sequence_checker(reference_id, reference_sequence, structure_id, structure_sequence, gapopen=10, gapextend=0.5,
                     allow_missing_on_termini=0, allow_mutants=False, allow_deletions=False,
                     allow_insertions=False, allow_unresolved=False,
                     write_output=False, outdir='', outfile='', force_rerun=False):
    """Report if a structure's sequence meets coverage checks to a reference sequence.

    First aligns a sequence from a chain of a PDB structure to "reference" sequence.
    Then creates a DataFrame of results and check for everything.

    Args:
        reference_id: ID of reference sequence
        reference_sequence: String representation of reference sequence
        structure_id: ID of sequence to be aligned
        structure_sequence: String representation of sequence to be aligned
        gapopen: Gap open penalty is the score taken away when a gap is created
        gapextend: Gap extension penalty is added to the standard gap penalty for each base or residue in the gap
        allow_missing_on_termini (float): Percentage of the length of the reference sequence which will be ignored
            when checking for modifications. Example: if 0.1, and reference sequence is 100 AA, then only residues
            10 to 90 will be checked for modifications.
        allow_mutants (bool): If mutations should be allowed or checked for
        allow_deletions (bool): If deletions should be allowed or checked for
        allow_insertions (bool): If insertions should be allowed or checked for
        allow_unresolved (bool): If unresolved residues should be allowed or checked for
        write_output (bool): Default False, set to True if you want the alignment file saved
        outdir (str, optional): Path to output directory. Default is the current directory.
        outfile (str, optional): Name of output file. If not set, is {id_a}_{id_b}_align.txt
        force_rerun (bool): Default False, set to True if you want to rerun the alignment if outfile exists.

    Returns:
        bool: If the structure's sequence meets the quality checks.

    """

    # Run the needle (global) alignment
    raw_alignment_results = ssbio.sequence.alignment.run_needle_alignment_on_str(id_a=reference_id,
                                                                                 seq_a=reference_sequence,
                                                                                 id_b=structure_id,
                                                                                 seq_b=structure_sequence,
                                                                                 gapopen=gapopen,
                                                                                 gapextend=gapextend,
                                                                                 write_output=write_output,
                                                                                 outdir=outdir,
                                                                                 outfile=outfile,
                                                                                 force_rerun=force_rerun)
    # Parse the alignment results
    summary_df = ssbio.sequence.alignment.get_alignment_summary_df(StringIO(raw_alignment_results))

    # Get cutoff stuff ready
    ref_seq_len = len(reference_sequence)
    # If any differences appear before start, they are ignored
    start = ref_seq_len - ref_seq_len * (1 - allow_missing_on_termini)
    # If any differences appear before end, they are ignored
    end = ref_seq_len - ref_seq_len * allow_missing_on_termini

    no_deletions_in_pdb = False
    no_insertions_in_pdb = False
    no_mutants_in_pdb = False
    no_unresolved_in_pdb = False

    # Check everything
    if not allow_deletions:
        # Get indices of the deletions
        deletion_indices = summary_df[summary_df['type'] == 'deletion'].index
        # If there are no deletions, that's great
        if len(deletion_indices) == 0:
            no_deletions_in_pdb = True
        else:
            # If the deletion appears before or after the cutoff, that's also great
            for deletion_index in deletion_indices:
                if summary_df.ix[deletion_index, 'id_a_stop'] < start or summary_df.ix[deletion_index, 'id_a_start'] > end:
                    no_deletions_in_pdb = True
                else:
                    no_deletions_in_pdb = False
                    log.debug('{} vs {}: PDB has deletions within structure'.format(reference_id, structure_id))
                    break
    else:
        no_deletions_in_pdb = True

    if not allow_insertions:
        # Get indices of the insertions
        insertion_indices = summary_df[summary_df['type'] == 'insertion'].index
        # If there are no insertions, that's great
        if len(insertion_indices) == 0:
            no_insertions_in_pdb = True
        else:
            # If the insertion appears before or after the cutoff, that's also great
            for insertion_index in insertion_indices:
                if summary_df.ix[insertion_index, 'id_a_stop'] < start or summary_df.ix[insertion_index, 'id_a_start'] > end:
                    no_insertions_in_pdb = True
                else:
                    no_insertions_in_pdb = False
                    log.debug('{} vs {}: PDB has insertions within structure'.format(reference_id, structure_id))
                    break
    else:
        no_insertions_in_pdb = True

    if not allow_mutants:
        # Get indices of the mutants
        mutant_indices = summary_df[summary_df['type'] == 'mutant'].index
        # If there are no mutants, that's great
        if len(mutant_indices) == 0:
            no_mutants_in_pdb = True
        else:
            # If the mutant appears before or after the cutoff, that's also great
            for mutant_index in mutant_indices:
                if summary_df.ix[mutant_index, 'id_a_stop'] < start or summary_df.ix[mutant_index, 'id_a_start'] > end:
                    no_mutants_in_pdb = True
                else:
                    no_mutants_in_pdb = False
                    log.debug('{} vs {}: PDB has mutants within structure'.format(reference_id, structure_id))
                    break
    else:
        no_mutants_in_pdb = True

    if not allow_unresolved:
        # Get indices of the unresolved residues
        unresolved_indices = summary_df[summary_df['type'] == 'unresolved'].index
        # If there are no unresolved, that's great
        if len(unresolved_indices) == 0:
            no_unresolved_in_pdb = True
        else:
            # If the unresolved residue appears before or after the cutoff, that's also great
            for unresolved_index in unresolved_indices:
                if summary_df.ix[unresolved_index, 'id_a_stop'] < start or summary_df.ix[unresolved_index, 'id_a_start'] > end:
                    no_unresolved_in_pdb = True
                else:
                    no_unresolved_in_pdb = False
                    log.debug('{} vs {}: PDB has unresolved within structure'.format(reference_id, structure_id))
                    break
    else:
        no_unresolved_in_pdb = True

    if no_deletions_in_pdb and no_insertions_in_pdb and no_mutants_in_pdb and no_unresolved_in_pdb:
        return True
    else:
        return False


def parse_procheck(quality_directory):
    """Parses all PROCHECK files in a directory and returns a Pandas DataFrame of the results

    Args:
        quality_directory: path to directory with PROCHECK output (.sum files)

    Returns:
        Pandas DataFrame: Summary of PROCHECK results

    """

    # TODO: save as dict instead, offer df as option
    # TODO: parse for one file instead

    procheck_summaries = glob.glob(os.path.join(quality_directory, '*.sum'))

    if len(procheck_summaries) == 0:
        return pd.DataFrame()

    all_procheck = {}
    for summ in procheck_summaries:
        structure_id = os.path.basename(summ).split('.sum')[0]
        procheck_dict = {}
        with open(summ) as f_in:
            lines = (line.rstrip() for line in f_in) # All lines including the blank ones
            lines = (line for line in lines if line) # Non-blank lines
            for line in lines:
                if len(line.split()) > 1:
                    if line.split()[1] == 'Ramachandran':
                        procheck_dict['procheck_rama_favored'] = percentage_to_float(line.split()[3])
                        procheck_dict['procheck_rama_allowed'] = percentage_to_float(line.split()[5])
                        procheck_dict['procheck_rama_allowed_plus'] = percentage_to_float(line.split()[7])
                        procheck_dict['procheck_rama_disallowed'] = percentage_to_float(line.split()[9])
                    if line.split()[1] == 'G-factors':
                        procheck_dict['procheck_gfac_dihedrals'] = line.split()[3]
                        procheck_dict['procheck_gfac_covalent'] = line.split()[5]
                        procheck_dict['procheck_gfac_overall'] = line.split()[7]
        all_procheck[structure_id] = procheck_dict

    DF_PROCHECK = pd.DataFrame.from_dict(all_procheck, orient='index')

    return DF_PROCHECK


def parse_psqs(psqs_results_file):
    """Parse a PSQS result file and returns a Pandas DataFrame of the results

    Args:
        psqs_results_file: Path to psqs results file

    Returns:
        Pandas DataFrame: Summary of PSQS results

    """

    # TODO: generalize column names for all results, save as dict instead

    psqs_results = pd.read_csv(psqs_results_file, sep='\t', header=None)
    psqs_results['pdb_file'] = psqs_results[0].apply(lambda x: str(x).strip('./').strip('.pdb'))
    psqs_results = psqs_results.rename(columns = {1:'psqs_local', 2:'psqs_burial', 3:'psqs_contact', 4:'psqs_total'}).drop(0, axis=1)
    psqs_results['u_pdb'] = psqs_results['pdb_file'].apply(lambda x: x.upper() if len(x)==4 else np.nan)
    psqs_results['i_entry_name'] = psqs_results['pdb_file'].apply(lambda x: x.split('_model1')[0] if len(x)>4 else np.nan)
    psqs_results = psqs_results[pd.notnull(psqs_results.psqs_total)]

    return psqs_results
