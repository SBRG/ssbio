import glob
import logging
import os

import numpy as np
import pandas as pd

import ssbio.protein.sequence.utils.alignment
from ssbio.utils import percentage_to_float

log = logging.getLogger(__name__)


def seq_to_struct_alignment_stats(reference_seq_aln, structure_seq_aln):
    """Get a report of an alignment from a sequence to a structure chain's sequence.

    Args:
        reference_seq_aln (str, Seq, SeqRecord): Reference sequence, alignment form
        structure_seq_aln (str, Seq, SeqRecord): Structure sequence, alignment form

    Returns:
        dict: Dictionary of information on mutations, insertions, sequence identity, etc.

    """
    if len(reference_seq_aln) != len(structure_seq_aln):
        raise ValueError('Sequence lengths not equal - was an alignment run?')

    reference_seq_aln = ssbio.protein.sequence.utils.cast_to_str(reference_seq_aln)
    structure_seq_aln = ssbio.protein.sequence.utils.cast_to_str(structure_seq_aln)

    infodict = {}

    # Percent identity to the reference sequence
    stats_percent_ident = ssbio.protein.sequence.utils.alignment.get_percent_identity(reference_seq_aln,
                                                                                      structure_seq_aln)
    infodict['percent_identity'] = stats_percent_ident

    # Other alignment results
    aln_df = ssbio.protein.sequence.utils.alignment.get_alignment_df(a_aln_seq=reference_seq_aln,
                                                                     b_aln_seq=structure_seq_aln)
    infodict['deletions'] = ssbio.protein.sequence.utils.alignment.get_deletions(aln_df)
    infodict['insertions'] = ssbio.protein.sequence.utils.alignment.get_insertions(aln_df)
    infodict['mutations'] = ssbio.protein.sequence.utils.alignment.get_mutations(aln_df)
    infodict['unresolved'] = ssbio.protein.sequence.utils.alignment.get_unresolved(aln_df)

    return infodict


def sequence_checker(reference_seq_aln, structure_seq_aln,
                     seq_ident_cutoff=0.5, allow_missing_on_termini=0.2,
                     allow_mutants=False, allow_deletions=False,
                     allow_insertions=False, allow_unresolved=False):
    """Report if a structure's sequence meets coverage checks to a reference sequence.

    First aligns a sequence from a chain of a PDB structure to "reference" sequence.
    Then creates a DataFrame of results to check for everything.

    Args:
        reference_seq_aln (str, Seq, SeqRecord): Reference sequence, alignment form
        structure_seq_aln (str, Seq, SeqRecord): Structure sequence, alignment form
        seq_ident_cutoff (float): Percent sequence identity cutoff, in decimal form
        allow_missing_on_termini (float): Percentage of the total length of the reference sequence which will be ignored
            when checking for modifications. Example: if 0.1, and reference sequence is 100 AA, then only residues
            5 to 95 will be checked for modifications.
        allow_mutants (bool): If mutations should be allowed or checked for
        allow_deletions (bool): If deletions should be allowed or checked for
        allow_insertions (bool): If insertions should be allowed or checked for
        allow_unresolved (bool): If unresolved residues should be allowed or checked for

    Returns:
        bool: If the structure's sequence meets the quality checks.

    """
    reference_seq_aln = ssbio.protein.sequence.utils.cast_to_str(reference_seq_aln)
    structure_seq_aln = ssbio.protein.sequence.utils.cast_to_str(structure_seq_aln)
    results = seq_to_struct_alignment_stats(reference_seq_aln=reference_seq_aln, structure_seq_aln=structure_seq_aln)

    # Check percent identity cutoff
    stats_percent_ident = results['percent_identity']
    log.debug('{}: percent identity'.format(stats_percent_ident))
    if stats_percent_ident < seq_ident_cutoff:
        log.debug('Alignment does not meet percent identity cutoff')
        return False
    else:
        log.debug('Alignment meets percent identity cutoff')

    # Get cutoff stuff ready
    ref_seq_len = len(reference_seq_aln.replace('-', ''))
    allow_missing_on_termini /= 2
    # If any differences appear before start, they are ignored
    start = ref_seq_len - (ref_seq_len * (1 - allow_missing_on_termini))
    # If any differences appear before end, they are ignored
    end = ref_seq_len - (ref_seq_len * allow_missing_on_termini)

    no_deletions_in_pdb = False
    no_insertions_in_pdb = False
    no_mutants_in_pdb = False
    no_unresolved_in_pdb = False

    # Check everything
    if not allow_deletions:
        # Get indices of the deletions
        deletions = results['deletions']
        # If there are no deletions, that's great
        if len(deletions) == 0:
            log.debug('No deletion regions')
            no_deletions_in_pdb = True
        else:
            log.debug('{} deletion region(s)'.format(len(deletions)))
            # If the deletion appears before or after the cutoff, that's also great
            for deletion in deletions:
                if deletion[0][1] < start or deletion[0][0] > end:
                    no_deletions_in_pdb = True
                    log.debug('Deletion region(s) are not within structure core')
                else:
                    no_deletions_in_pdb = False
                    log.debug('Deletions within structure')
                    log.debug('{} > {} or {} < {}'.format(deletion[0][1], start, deletion[0][0], end))
                    break
    else:
        no_deletions_in_pdb = True

    if not allow_insertions:
        # Get indices of the insertions
        insertions = results['insertions']
        # If there are no insertions, that's great
        if len(insertions) == 0:
            log.debug('No insertion regions')
            no_insertions_in_pdb = True
        else:
            log.debug('{} insertion region(s)'.format(len(insertions)))
            # If the insertion appears before or after the cutoff, that's also great
            for insertion in insertions:
                if insertion[0][1] < start or insertion[0][0] > end:
                    no_insertions_in_pdb = True
                    log.debug('Insertion region(s) are not within structure core')
                else:
                    no_insertions_in_pdb = False
                    log.debug('Insertion regions within structure')
                    break
    else:
        no_insertions_in_pdb = True

    if not allow_mutants:
        # Get indices of the mutants
        mutations_full = results['mutations']
        mutations = [x[1] for x in mutations_full]
        # If there are no mutants, that's great
        if len(mutations) == 0:
            log.debug('No point mutations')
            no_mutants_in_pdb = True
        else:
            log.debug('{} point mutation(s)'.format(len(mutations)))
            # If the mutant appears before or after the cutoff, that's also great
            for mutation in mutations:
                if mutation < start or mutation > end:
                    no_mutants_in_pdb = True
                    log.debug('Mutation region(s) are not within structure core')
                else:
                    no_mutants_in_pdb = False
                    log.debug('Mutantion regions within structure')
                    break
    else:
        no_mutants_in_pdb = True

    if not allow_unresolved:
        # Get indices of the unresolved residues
        unresolved = results['unresolved']
        # If there are no unresolved, that's great
        if len(unresolved) == 0:
            log.debug('No unresolved mutations')
            no_unresolved_in_pdb = True
        else:
            log.debug('{} unresolved residue(s)'.format(len(unresolved)))
            # If the unresolved residue appears before or after the cutoff, that's also great
            for unr in unresolved:
                if unr < start or unr > end:
                    no_unresolved_in_pdb = True
                    log.debug('Unresolved region(s) are not within structure core')
                else:
                    no_unresolved_in_pdb = False
                    log.debug('Unresolved residues within structure')
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
