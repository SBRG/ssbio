import re
import os.path as op
import cobra.io.dict
from cobra.io import load_matlab_model
from cobra.io import read_sbml_model
from cobra.io import load_json_model
from cobra.core import Model
from cobra.core import Gene
from cobra.core import Reaction
from cobra.core import DictList
from ssbio.core.genepro import GenePro
import numpy as np
import networkx as nx


class ModelPro(Model):
    def __new__(cls, *args, **kwargs):
        """Casting Gene objects into GenePro objects
            This replaces any new instance of "X" with "X"-Pros. Even when you load a model later
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

    def store_network_representations(self, exclude_mets=None):
        # Remove highly connected metabolites
        excluded_mets_tmp = ['atp', 'adp', 'amp', 'pi', 'ppi', 'co2',
                             'h', 'nad', 'nadp', 'coa', 'glu__L', 'glu__D',
                             'nadph', 'nadphx__R', 'nadphx__S', 'h2o', 'nh4', 'nadh']
        excluded_mets = []
        for x in excluded_mets_tmp:
            excluded_mets.append(x + '_c')
            excluded_mets.append(x + '_e')
            excluded_mets.append(x + '_p')

        if exclude_mets:
            assert isinstance(exclude_mets, list)
            excluded_mets.extend(exclude_mets)

        # Create gene to metabolite (product or reactant) dictionary
        g_to_met_dict = {}
        for g in self.genes:
            g_rxns = [r for r in list(g.reactions)]
            g_mets = []
            for r in g_rxns:
                for m in r.metabolites:
                    if m.id not in excluded_mets:
                        g_mets.append(m.id)
            g_to_met_dict[g.id] = g_mets

        # Genes as nodes, metabolites as edges
        GG = nx.Graph()
        for g1, g1_mets in g_to_met_dict.items():
            for g2, g2_mets in g_to_met_dict.items():
                if g1 != g2:
                    incommon = list(set(g1_mets).intersection(g2_mets))
                    if len(incommon) > 0:
                        # Just create edge if there's anything in common
                        GG.add_edge(g1, g2, object=incommon)
        self.gene_met_network = GG

        # Metabolites as nodes, reactions as edges
        GM = nx.DiGraph()
        for r in self.reactions:
            # Can carry forward flux?
            if r.upper_bound > 0:
                for reactant in r.reactants:
                    for product in r.products:
                        if reactant.id in excluded_mets or product.id in excluded_mets:
                            continue
                        GM.add_edge(reactant.id, product.id, object=r.id)
            # Can carry reverse flux?
            if r.lower_bound > 0:
                for reactant in r.reactants:
                    for product in r.products:
                        if reactant.id in excluded_mets or product.id in excluded_mets:
                            continue
                        GM.add_edge(product.id, reactant.id, object=r.id)
        self.met_rxn_network = GM

    def network_distance_between_genes(self, gene1_id, gene2_id, shortest=True):
        all_shortest_paths = []
        for i, path in enumerate(nx.all_shortest_paths(self.gene_met_network, gene1_id, gene2_id)):
            edges_in_path = zip(path[0:], path[1:])
            full_info = [(u, v, self.gene_met_network[u][v]['object']) for (u, v) in edges_in_path]
            if shortest:
                return full_info
            all_shortest_paths.append(full_info)
        return all_shortest_paths

    def network_distance_between_mets(self, met1_id, met2_id, shortest=True):
        all_shortest_paths = []
        for i, path in enumerate(nx.all_shortest_paths(self.met_rxn_network, met1_id, met2_id)):
            edges_in_path = zip(path[0:], path[1:])
            full_info = [(u, v, self.met_rxn_network[u][v]['object']) for (u, v) in edges_in_path]
            if shortest:
                return full_info
            all_shortest_paths.append(full_info)
        return all_shortest_paths


def model_loader(gem_file_path, gem_file_type):
    """Consolidated function to load a GEM using COBRApy. Specify the file type being loaded.

    Args:
        gem_file_path (str): Path to model file
        gem_file_type (str): GEM model type - ``sbml`` (or ``xml``), ``mat``, or ``json`` format

    Returns:
        COBRApy Model object.

    """

    if gem_file_type.lower() == 'xml' or gem_file_type.lower() == 'sbml':
        model = read_sbml_model(gem_file_path)
    elif gem_file_type.lower() == 'mat':
        model = load_matlab_model(gem_file_path)
    elif gem_file_type.lower() == 'json':
        model = load_json_model(gem_file_path)
    else:
        raise ValueError('File type must be "sbml", "xml", "mat", or "json".')

    return model


def is_spontaneous(gene, custom_id=None):
    """Input a COBRApy Gene object and check if the ID matches a spontaneous ID regex.
    
    Args:
        gene (Gene): COBRApy Gene
        custom_id (str): Optional custom spontaneous ID if it does not match the regular expression ``[Ss](_|)0001``

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
        genes (DictList): Genes DictList
        custom_spont_id (str): Optional custom spontaneous ID if it does not match the regular expression ``[Ss](_|)0001``

    Returns:
        DictList: genes excluding ones that are spontaneous

    """
    new_genes = DictList()
    for gene in genes:
        if not is_spontaneous(gene, custom_id=custom_spont_id):
            new_genes.append(gene)

    return new_genes


def true_num_genes(model, custom_spont_id=None):
    """Return the number of genes in a model ignoring spontaneously labeled genes.

    Args:
        model (Model):
        custom_spont_id (str): Optional custom spontaneous ID if it does not match the regular expression ``[Ss](_|)0001``

    Returns:
        int: Number of genes excluding spontaneous genes

    """
    true_num = 0
    for gene in model.genes:
        if not is_spontaneous(gene, custom_id=custom_spont_id):
            true_num += 1
    return true_num


def true_num_reactions(model, custom_spont_id=None):
    """Return the number of reactions associated with a gene.

    Args:
        model (Model):
        custom_spont_id (str): Optional custom spontaneous ID if it does not match the regular expression ``[Ss](_|)0001``

    Returns:
        int: Number of reactions associated with a gene

    """
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