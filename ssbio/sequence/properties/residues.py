from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.Polypeptide import one_to_three
import ssbio.utils
import ssbio.sequence.utils
import logging
log = logging.getLogger(__name__)


# TODO: function to predict secondary structure (PSIPRED)
# TODO: function to predict solvent accessibility/resdepth just based on sequence (SCRATCH, RDPred)
# See more: https://www.researchgate.net/publication/235633540_Recent_Advances_in_Predicting_Functional_Impact_of_Single_Amino_Acid_Polymorphisms_A_Review_of_Useful_Features_Computational_Methods_and_Available_Tools

_aa_property_dict_one = {
    'Aliphatic': ['A', 'I', 'L', 'V'],
    'Aromatic' : ['F', 'H', 'W', 'Y'],
    'Non-polar': ['A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y'],
    'Polar'    : ['D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T'],
    'Charged'  : ['D', 'E', 'H', 'K', 'R'],
    'Basic'    : ['H', 'K', 'R'],
    'Acidic'   : ['D', 'E']}
# 'Tiny': ['A','C','G','S','T']
# 'Small': ['A','C','D','G','N','P','S','T','V']

_aa_property_dict_three = {k: [one_to_three(x) for x in v] for k, v in _aa_property_dict_one.items()}


def biopython_protein_analysis(inseq):
    """Utiize Biopython's ProteinAnalysis module to return general sequence properties of an amino acid string.

    For full definitions see: http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam.ProteinAnalysis-class.html

    Args:
        inseq: Amino acid sequence

    Returns:
        dict: Dictionary of sequence properties. Some definitions include:
        instability_index: Any value above 40 means the protein is unstable (has a short half life).
        secondary_structure_fraction: Percentage of protein in helix, turn or sheet

    TODO:
        Finish definitions of dictionary

    """

    inseq = ssbio.sequence.utils.cast_to_str(inseq)

    analysed_seq = ProteinAnalysis(inseq)

    info_dict = {}
    # info_dict['amino_acids_content'] = analysed_seq.count_amino_acids()
    # info_dict['amino_acids_percent'] = analysed_seq.get_amino_acids_percent()
    # info_dict['length'] = analysed_seq.length
    info_dict['monoisotopic'] = analysed_seq.monoisotopic
    info_dict['molecular_weight'] = analysed_seq.molecular_weight()
    info_dict['aromaticity'] = analysed_seq.aromaticity()
    info_dict['instability_index'] = analysed_seq.instability_index()
    # TODO: What is flexibility?
    # info_dict['flexibility'] = analysed_seq.flexibility()
    info_dict['isoelectric_point'] = analysed_seq.isoelectric_point()

    # Separated secondary_structure_fraction into each definition
    # info_dict['secondary_structure_fraction'] = analysed_seq.secondary_structure_fraction()
    info_dict['percent_helix_naive'] = analysed_seq.secondary_structure_fraction()[0]
    info_dict['percent_turn_naive'] = analysed_seq.secondary_structure_fraction()[1]
    info_dict['percent_sheet_naive'] = analysed_seq.secondary_structure_fraction()[2]

    return info_dict


def emboss_pepstats_on_fasta(infile, outfile='', outdir='', outext='.pepstats', force_rerun=False):
    """Run EMBOSS pepstats on a FASTA file.

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
    outfile = ssbio.utils.outfile_maker(inname=infile, outname=outfile, outdir=outdir, outext=outext)

    # Run pepstats
    program = 'pepstats'
    pepstats_args = '-sequence="{}" -outfile="{}"'.format(infile, outfile)
    cmd_string = '{} {}'.format(program, pepstats_args)
    ssbio.utils.command_runner(cmd_string, force_rerun_flag=force_rerun, outfile=outfile, silent=True)

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


def residue_biochemical_definition(res):
    # TODO: docstring
    resprop = []
    for k, v in _aa_property_dict_one.items():
        if res in v:
            resprop.append(k)

    return resprop


def grantham_score(ref_aa, mut_aa):
    """https://github.com/ashutoshkpandey/Annotation/blob/master/Grantham_score_calculator.py"""
    grantham = {
        'S': {'R': 110, 'L': 145, 'P': 74, 'T': 58, 'A': 99, 'V': 124, 'G': 56, 'I': 142, 'F': 155, 'Y': 144, 'C': 112,
              'H': 89, 'Q': 68, 'N': 46, 'K': 121, 'D': 65, 'E': 80, 'M': 135, 'W': 177},
        'R': {'R': 0, 'L': 102, 'P': 103, 'T': 71, 'A': 112, 'V': 96, 'G': 125, 'I': 97, 'F': 97, 'Y': 77, 'C': 180,
              'H': 29, 'Q': 43, 'N': 86, 'K': 26, 'D': 96, 'E': 54, 'M': 91, 'W': 101, 'S': 0},
        'L': {'R': 0, 'L': 0, 'P': 98, 'T': 92, 'A': 96, 'V': 32, 'G': 138, 'I': 5, 'F': 22, 'Y': 36, 'C': 198, 'H': 99,
              'Q': 113, 'N': 153, 'K': 107, 'D': 172, 'E': 138, 'M': 15, 'W': 61, 'S': 0},
        'P': {'R': 0, 'L': 0, 'P': 0, 'T': 38, 'A': 27, 'V': 68, 'G': 42, 'I': 95, 'F': 114, 'Y': 110, 'C': 169,
              'H': 77, 'Q': 76, 'N': 91, 'K': 103, 'D': 108, 'E': 93, 'M': 87, 'W': 147, 'S': 0},
        'T': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 58, 'V': 69, 'G': 59, 'I': 89, 'F': 103, 'Y': 92, 'C': 149, 'H': 47,
              'Q': 42, 'N': 65, 'K': 78, 'D': 85, 'E': 65, 'M': 81, 'W': 128, 'S': 0},
        'A': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 64, 'G': 60, 'I': 94, 'F': 113, 'Y': 112, 'C': 195, 'H': 86,
              'Q': 91, 'N': 111, 'K': 106, 'D': 126, 'E': 107, 'M': 84, 'W': 148, 'S': 0},
        'V': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 109, 'I': 29, 'F': 50, 'Y': 55, 'C': 192, 'H': 84,
              'Q': 96, 'N': 133, 'K': 97, 'D': 152, 'E': 121, 'M': 21, 'W': 88, 'S': 0},
        'G': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 0, 'I': 135, 'F': 153, 'Y': 147, 'C': 159, 'H': 98,
              'Q': 87, 'N': 80, 'K': 127, 'D': 94, 'E': 98, 'M': 127, 'W': 184, 'S': 0},
        'I': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 0, 'I': 0, 'F': 21, 'Y': 33, 'C': 198, 'H': 94,
              'Q': 109, 'N': 149, 'K': 102, 'D': 168, 'E': 134, 'M': 10, 'W': 61, 'S': 0},
        'F': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 0, 'I': 0, 'F': 0, 'Y': 22, 'C': 205, 'H': 100,
              'Q': 116, 'N': 158, 'K': 102, 'D': 177, 'E': 140, 'M': 28, 'W': 40, 'S': 0},
        'Y': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 0, 'I': 0, 'F': 0, 'Y': 0, 'C': 194, 'H': 83,
              'Q': 99, 'N': 143, 'K': 85, 'D': 160, 'E': 122, 'M': 36, 'W': 37, 'S': 0},
        'C': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 0, 'I': 0, 'F': 0, 'Y': 0, 'C': 0, 'H': 174,
              'Q': 154, 'N': 139, 'K': 202, 'D': 154, 'E': 170, 'M': 196, 'W': 215, 'S': 0},
        'H': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 0, 'I': 0, 'F': 0, 'Y': 0, 'C': 0, 'H': 0, 'Q': 24,
              'N': 68, 'K': 32, 'D': 81, 'E': 40, 'M': 87, 'W': 115, 'S': 0},
        'Q': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 0, 'I': 0, 'F': 0, 'Y': 0, 'C': 0, 'H': 0, 'Q': 0,
              'N': 46, 'K': 53, 'D': 61, 'E': 29, 'M': 101, 'W': 130, 'S': 0},
        'N': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 0, 'I': 0, 'F': 0, 'Y': 0, 'C': 0, 'H': 0, 'Q': 0,
              'N': 0, 'K': 94, 'D': 23, 'E': 42, 'M': 142, 'W': 174, 'S': 0},
        'K': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 0, 'I': 0, 'F': 0, 'Y': 0, 'C': 0, 'H': 0, 'Q': 0,
              'N': 0, 'K': 0, 'D': 101, 'E': 56, 'M': 95, 'W': 110, 'S': 0},
        'D': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 0, 'I': 0, 'F': 0, 'Y': 0, 'C': 0, 'H': 0, 'Q': 0,
              'N': 0, 'K': 0, 'D': 0, 'E': 45, 'M': 160, 'W': 181, 'S': 0},
        'E': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 0, 'I': 0, 'F': 0, 'Y': 0, 'C': 0, 'H': 0, 'Q': 0,
              'N': 0, 'K': 0, 'D': 0, 'E': 0, 'M': 126, 'W': 152, 'S': 0},
        'M': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 0, 'I': 0, 'F': 0, 'Y': 0, 'C': 0, 'H': 0, 'Q': 0,
              'N': 0, 'K': 0, 'D': 0, 'E': 0, 'M': 0, 'W': 67, 'S': 0},
        'W': {'R': 0, 'L': 0, 'P': 0, 'T': 0, 'A': 0, 'V': 0, 'G': 0, 'I': 0, 'F': 0, 'Y': 0, 'C': 0, 'H': 0, 'Q': 0,
              'N': 0, 'K': 0, 'D': 0, 'E': 0, 'M': 0, 'W': 0, 'S': 0}}

    score = 0

    if ref_aa not in grantham or mut_aa not in grantham:
        log.error('{} to {}: a residue is not in the Grantham matrix'.format(ref_aa, mut_aa))
        return score, 'Unknown'

    if ref_aa == mut_aa:
        return score, 'Conservative'
    else:
        if int(grantham[ref_aa][mut_aa]) != 0:
            score = score + int(grantham[ref_aa][mut_aa])
        else:
            score = score + int(grantham[mut_aa][ref_aa])

    if score > 150:
        return score, "Radical"
    elif score <= 150 and score > 100:
        return score, "Moderately Radical"
    elif score <= 100 and score > 50:
        return score, "Moderately Conservative"
    else:
        return score, "Conservative"

