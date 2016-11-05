from Bio.SeqUtils.ProtParam import ProteinAnalysis
import ssbio.utils
import subprocess
import logging
log = logging.getLogger(__name__)


def sequence_properties(seq_str):
    """Utiize Biopython's ProteinAnalysis module to return general sequence properties of an amino acid string.

    Args:
        seq_str: String representation of a amino acid sequence

    Returns:
        dict: Dictionary of sequence properties. Some definitions include:
        instability_index: Any value above 40 means the protein is unstable (has a short half life).
        secondary_structure_fraction: Percentage of protein in helix, turn or sheet

    TODO:
        Finish definitions of dictionary

    """

    analysed_seq = ProteinAnalysis(seq_str)

    info_dict = {}
    info_dict['amino_acids_content'] = analysed_seq.count_amino_acids()
    info_dict['amino_acids_percent'] = analysed_seq.get_amino_acids_percent()
    info_dict['length'] = analysed_seq.length
    info_dict['monoisotopic'] = analysed_seq.monoisotopic
    info_dict['molecular_weight'] = analysed_seq.molecular_weight()
    info_dict['aromaticity'] = analysed_seq.aromaticity()
    info_dict['instability_index'] = analysed_seq.instability_index()
    info_dict['flexibility'] = analysed_seq.flexibility()
    info_dict['isoelectric_point'] = analysed_seq.isoelectric_point()
    info_dict['secondary_structure_fraction'] = analysed_seq.secondary_structure_fraction()

    return info_dict


def emboss_pepstats_on_fasta(infile, outfile='', outdir='', outext='.pepstats', force_rerun=False):
    """Run EMBOSS pepstats on a sequence string, or FASTA file.

    Args:
        infile: Path to FASTA file
        outfile: Name of output file without extension
        outdir: Path to output directory
        outext: Extension of results file, default is ".pepstats"
        force_rerun: Flag to rerun pepstats

    Returns:
        str: Path to output file.

    """

    # Create the output file name
    outfile = ssbio.utils.outfile_name_maker(inname=infile, outfile=outfile, outdir=outdir, outext=outext)

    # Run pepstats
    pepstats_args = '-sequence="{}" -outfile="{}"'.format(infile, outfile)
    ssbio.utils.command_runner(program='pepstats', args=pepstats_args, force_rerun_flag=force_rerun, outfile=outfile)

    return outfile


def emboss_pepstats_on_str(instring, outfile, outdir='', outext='.pepstats', force_rerun=False):
    """Run EMBOSS pepstats on a sequence string, or FASTA file.

    Args:
        instring: Sequence string
        outfile: Name of output file without extension
        outdir: Path to output directory
        outext: Extension of results file, default is ".pepstats"
        force_rerun: Flag to rerun pepstats

    Returns:
        str: Path to output file.

    """
    # Create the output file name
    outfile = ssbio.utils.outfile_name_maker(inname='seq_str', outfile=outfile, outdir=outdir, outext=outext)

    # Run pepstats
    pepstats_args = '-sequence=asis::{} -outfile="{}"'.format(instring, outfile)
    ssbio.utils.command_runner(program='pepstats', args=pepstats_args, force_rerun_flag=force_rerun, outfile=outfile)

    return outfile


def emboss_pepstats_parser(infile):
    """Get dictionary of pepstats results.

    Args:
        infile: Path to pepstats outfile

    Returns:
        dict: Parsed information from pepstats

    TODO:
        Only currently parsing the bottom of the file for percentages of properties.

    """
    with open(infile) as f:
        lines = f.read().split('\n')

    info_dict = {}

    for l in lines[38:47]:
        info = l.split('\t')
        cleaninfo = list(filter(lambda x: x != '', info))
        prop = cleaninfo[0]
        num = cleaninfo[2]
        percent = float(cleaninfo[-1]) / float(100)

        info_dict['percent_' + prop.lower()] = percent

    return info_dict

#
# AAdict = {
#
#
# 'LYS': 'positive',
# 'ARG': 'positive',
# 'HIS': 'positive',
#
# 'ASP': 'negative',
# 'GLU': 'negative',
#
# 'LEU': 'nonpolar',
# 'TRP': 'nonpolar',
# 'VAL': 'nonpolar',
# 'PHE': 'nonpolar',
# 'PRO': 'nonpolar',
# 'ILE': 'nonpolar',
# 'GLY': 'nonpolar',
# 'ALA': 'nonpolar',
# 'MET': 'nonpolar',
#
# 'ASN': 'polar',
# 'THR': 'polar',
# 'TYR': 'polar',
# 'MSE': 'polar',
# 'SEC': 'polar',
# 'SER': 'polar',
# 'GLN': 'polar',
# 'CYS': 'polar',
# }
#
#
# def residue_props(seq_str):
#     """Return a dictionary of residue properties indicating the percentage of the respective property for a sequence.
#
#     Properties are: Polar, nonpolar, negative, positive.
#
#     Args:
#         pdb_file: PDB or MMCIF structure file
#
#     Returns:
#         dict: Dictonary of percentage (float) of properties
#     """
#
#     props = {}
#     polar = 0
#     nonpolar = 0
#     positive = 0
#     negative = 0
#     total = 0
#
#     for j in seq_str:
#         if j.resname in AAdict:
#             if AAdict[j.resname] == 'nonpolar':
#                 nonpolar = nonpolar + 1
#             elif AAdict[j.resname] == 'polar':
#                 polar = polar + 1
#             elif AAdict[j.resname] == 'positive':
#                 positive = positive + 1
#             elif AAdict[j.resname] == 'negative':
#                 negative = negative + 1
#             total = total + 1
#
#     props['ssb_per_NP'] = float(nonpolar) / float(total)
#     props['ssb_per_P'] = float(polar) / float(total)
#     props['ssb_per_pos'] = float(positive) / float(total)
#     props['ssb_per_neg'] = float(negative) / float(total)
#
#     return props