import re
import os.path as op
import cobra.io.dict
from cobra.io import load_matlab_model
from cobra.io import read_sbml_model
from cobra.io import load_json_model
from cobra.core import Model
from cobra.core import Gene
from cobra.core import DictList
from ssbio.core.genepro import GenePro


class ModelPro(Model):
    def __new__(cls, *args, **kwargs):
        """Casting Gene objects into GenePro objects
            This replaces any new instance of Genes with GenePros. Even when you load a model later
            using COBRApy methods. Use with caution!
        See http://stackoverflow.com/questions/3464061/cast-base-class-to-derived-class-python-or-more-pythonic-way-of-extending-class

        Returns:
            GenePro: a Gene object with a .protein attribute
        """
        if cls == Gene:
            return object.__new__(GenePro)
        return object.__new__(cls)

    Gene.__new__ = staticmethod(__new__)

    def __json_encode__(self):
        """convert the model to a dict"""
        obj = dict(
                reactions=[cobra.io.dict.reaction_to_dict(reaction) for reaction in self.reactions],
                metabolites=[
                    cobra.io.dict.metabolite_to_dict(metabolite) for metabolite in self.metabolites
                    ],
                genes = self.genes,
                id=self.id,
        )
        cobra.io.dict._update_optional(self, obj, cobra.io.dict._OPTIONAL_MODEL_ATTRIBUTES, cobra.io.dict._ORDERED_OPTIONAL_MODEL_KEYS)
        # add in the JSON version
        obj["version"] = 1
        return obj

    def __json_decode__(self, **attrs):
        """build a model from a dict"""
        Model.__init__(self)
        if 'reactions' not in attrs:
            raise Exception('JSON object has no reactions attribute. Cannot load.')
        self.add_metabolites(
                [cobra.io.dict.metabolite_from_dict(metabolite) for metabolite in attrs['metabolites']]
        )
        self.genes = DictList(attrs['genes'])
        self.add_reactions(
                [cobra.io.dict.reaction_from_dict(reaction, self) for reaction in attrs['reactions']]
        )
        for k, v in attrs.items():
            if k in {'id', 'name', 'notes', 'compartments', 'annotation'}:
                setattr(self, k, v)


def model_loader(gem_file_path, gem_file_type, pdb_file_type='cif'):
    """Consolidated function to load a GEM using COBRApy. Specify the file type being loaded.

    Args:
        gem_file_path (str): Path to model file
        gem_file_type (str): GEM model type - 'sbml' (or 'xml'), 'mat', or 'json' format

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

    modelpro = ModelPro(model)
    for g in modelpro.genes:
        g.protein.pdb_file_type = pdb_file_type

    return modelpro


def is_spontaneous(gene, custom_id=None):
    """Input a COBRApy Gene object and check if the ID matches a spontaneous ID regex.
    
    Args:
        gene (Gene): COBRApy Gene
        custom_id (str): Custom spontaneous ID if not matching regex

    Returns:
        bool: If gene ID matches spontaneous ID

    """
    spont = re.compile("[Ss](_|)0001")
    if spont.match(gene.id):
        return True
    elif gene.id == custom_id:
        return True
    else:
        return False


def filter_out_spontaneous_genes(genes, custom_spont_id=None):
    """Return the DictList of genes that are not spontaneous in a model.

    Args:
        model: COBRApy Model object

    Returns:
        DictList: genes excluding ones that are spontaneous

    """
    new_genes = DictList()
    for gene in genes:
        if not is_spontaneous(gene, custom_id=custom_spont_id):
            new_genes.append(gene)

    return new_genes

def true_num_genes(model, custom_spont_id=None):
    """Return the number of genes in a model ignoring spontaneously labeled genes

    :param model: COBRApy Model object
    :return int: Number of genes excluding spontaneous genes
    """
    true_num = 0
    for gene in model.genes:
        if not is_spontaneous(gene, custom_id=custom_spont_id):
            true_num += 1
    return true_num


def true_num_reactions(model, custom_spont_id=None):
    true_num = 0
    for rxn in model.reactions:
        if len(rxn.genes) == 0:
            continue
        if len(rxn.genes) == 1 and is_spontaneous(list(rxn.genes)[0], custom_id=custom_spont_id):
            continue
        else:
            true_num += 1
    return true_num


def adj_num_reactions(model, missing_genes, custom_spont_id=None):
    adj_num = 0
    for rxn in model.reactions:
        if len(rxn.genes) == 0:
            continue
        if len(rxn.genes) == 1 and (is_spontaneous(list(rxn.genes)[0], custom_id=custom_spont_id) or list(rxn.genes)[0] in missing_genes):
            continue
        else:
            adj_num += 1
    return adj_num