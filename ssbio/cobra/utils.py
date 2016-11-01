import re
import os.path as op
from cobra.io import load_matlab_model
from cobra.io import read_sbml_model
from cobra.io import load_json_model


def model_loader(gem_file_path, gem_file_type):
    """Consolidated function to load a GEM using COBRApy. Specify the file type being loaded.

    Args:
        gem_file_path (str): Path to model file
        gem_file_type (str): if your model is in "sbml" (or "xml"), "mat", or 'json" format

    Returns:
        COBRApy Model object.

    """

    extension = op.splitext(gem_file_path)[1]

    if gem_file_type.replace('.', '').lower() == 'xml' or gem_file_type.replace('.', '').lower() == 'sbml':
        model = read_sbml_model(gem_file_path)
    elif extension.replace('.', '').lower() == 'mat':
        model = load_matlab_model(gem_file_path)
    elif extension.replace('.', '').lower() == 'json':
        model = load_json_model(gem_file_path)
    else:
        raise ValueError('File type must be in SBML, MATLAB, or JSON formats.')

    return model


def is_spontaneous(gene):
    """Input a COBRApy Gene object and check if the ID matches a spontaneous ID regex.

    :param gene:
    :return:
    """
    spont = re.compile("[Ss](_|)0001")
    if spont.match(gene.id):
        return True
    else:
        return False


def true_num_genes(model):
    """Return the number of genes in a model ignoring spontaneously labeled genes

    :param model: COBRApy Model object
    :return int: Number of genes excluding spontaneous genes
    """
    true_num = 0
    for gene in model.genes:
        if not is_spontaneous(gene):
            true_num += 1
    return true_num


def true_num_reactions(model):
    true_num = 0
    for rxn in model.reactions:
        if len(rxn.genes) == 0:
            continue
        if len(rxn.genes) == 1 and is_spontaneous(list(rxn.genes)[0]):
            continue
        else:
            true_num += 1
    return true_num


def adj_num_reactions(model, missing_genes):
    adj_num = 0
    for rxn in model.reactions:
        if len(rxn.genes) == 0:
            continue
        if len(rxn.genes) == 1 and (is_spontaneous(list(rxn.genes)[0]) or list(rxn.genes)[0] in missing_genes):
            continue
        else:
            adj_num += 1
    return adj_num


# def gene_name_to_id(full_model_xml):
#     '''
#     Input: the .xml filename of the 'full' model (in my case Recon2) which has the gene names and
#     associated gene IDs listed as species
#     Output: a dictionary with the following structure - {GeneName: GeneID, ...}
#     '''
#
#     full_model_sbml = reader.readSBML(full_model_xml)
#     m = full_model_sbml.getModel()
#
#     gene_dict = defaultdict(dict)
#
#     # 'species' includes the genes
#     for i in m.getListOfSpecies():
#         entrez = i.getId()
#
#         # some start with _NM
#         if entrez.startswith('_NM'):
#             name = i.getName()
#             annotation = i.getAnnotation()
#             gene_dict[name] = entrez
#
#         # all other genes start with _
#         elif entrez.startswith('_'):
#             name = i.getName()
#             annotation = i.getAnnotation()
#             idz = entrez.split('_')
#             idz.pop(0)
#
#             # accounting for genes that are connected
#             newidz = idz[0::3]
#             namez = name.split(':')
#             for x in range(len(newidz)):
#                 if namez[x] == 'null':
#                     continue
#                 gene_dict[namez[x]] = newidz[x]
#
#     return gene_dict
#
# def full_to_reduced(full_dict, reduced_list_of_genes):
#     '''
#     Input:
#         full_dict: the dictionary returned by function nameToId
#         reduced_list_of_genes: the list of genes of your reduced model (in my case iAB-RBC)
#     Output:
#         gene_mapping: the gene mapping for your reduced model (a dictionary of dictionaries with form {GeneName: {'ENTREZ': GeneID, 'UNIPROT': UniprotAcc#}, ...})
#         missing_in_full: genes that were not in the full model but in the reduced - stuff that needs to be looked at manually
#     '''
#     missing_in_full = []
#     gene_mapping = {}
#
#     for gene in reduced_list_of_genes:
#         genename = gene.name.strip('|&').split('.')[0].upper()
#         if genename in full_dict:
#             gene_mapping[genename] = full_dict[genename]
#         else:
#             missing_in_full.append(genename)
#
#     return gene_mapping, list(set(filter(None, missing_in_full)))