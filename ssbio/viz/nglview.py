import logging
import ssbio.utils
import seaborn as sns
log = logging.getLogger(__name__)


def add_residues_highlight_to_nglview(view, structure_resnums, chain, res_color='red'):
    """Add a residue number or numbers to an NGLWidget view object.

    Args:
        view (NGLWidget): NGLWidget view object
        structure_resnums (int, list): Residue number(s) to highlight, structure numbering
        chain (str, list): Chain ID or IDs of which residues are a part of. If not provided, all chains in the
            mapped_chains attribute will be used. If that is also empty, and exception is raised.
        res_color (str): Color to highlight residues with

    """
    chain = ssbio.utils.force_list(chain)

    if isinstance(structure_resnums, list):
        structure_resnums = list(set(structure_resnums))
    elif isinstance(structure_resnums, int):
        structure_resnums = ssbio.utils.force_list(structure_resnums)
    else:
        raise ValueError('Input must either be a residue number of a list of residue numbers')

    to_show_chains = '( '
    for c in chain:
        to_show_chains += ':{} or'.format(c)
    to_show_chains = to_show_chains.strip(' or ')
    to_show_chains += ' )'

    to_show_res = '( '
    for m in structure_resnums:
        to_show_res += '{} or '.format(m)
    to_show_res = to_show_res.strip(' or ')
    to_show_res += ' )'

    log.info('Selection: {} and not hydrogen and {}'.format(to_show_chains, to_show_res))

    view.add_ball_and_stick(selection='{} and not hydrogen and {}'.format(to_show_chains, to_show_res), color=res_color)


def add_scaled_residues_highlight_to_nglview(view, structure_resnums, chain_id, color='red',
                                             unique_colors=False, opacity_range=(0.5,1), scale_range=(.7, 10)):
    """Add a list of residue numbers (which may contain repeating residues) to a view, or add a dictionary of
        residue numbers to counts. Size and opacity of added residues are scaled by counts.

    Args:
        view (NGLWidget): NGLWidget view object
        structure_resnums (int, list, dict): Residue number(s) to highlight, or a dictionary of residue number to
            frequency count
        chain_id (str, list): Chain ID or IDs of which residues are a part of.
        color (str): Color to highlight residues with
        unique_colors (bool): If each mutation should be colored uniquely (will override color argument)
        opacity_range (tuple): Min/max opacity values (residues that have higher frequency counts will be opaque)
        scale_range (tuple): Min/max size values (residues that have higher frequency counts will be bigger)

    """
    # TODO: likely to move these functions to a separate nglview/utils folder since they are not coupled to the structure
    # TODO: add color by letter_annotations!
    chain_id = ssbio.utils.force_list(chain_id)

    if isinstance(structure_resnums, dict):
        opacity_dict = ssbio.utils.scale_calculator(opacity_range[0], structure_resnums, rescale=opacity_range)
        scale_dict = ssbio.utils.scale_calculator(scale_range[0], structure_resnums, rescale=scale_range)
    else:
        opacity_dict = {x: max(opacity_range) for x in ssbio.utils.force_list(structure_resnums)}
        scale_dict = {x: max(scale_range) for x in ssbio.utils.force_list(structure_resnums)}

    if isinstance(structure_resnums, list):
        structure_resnums = list(set(structure_resnums))
    elif isinstance(structure_resnums, dict):
        structure_resnums = list(structure_resnums.keys())
    elif isinstance(structure_resnums, int):
        structure_resnums = ssbio.utils.force_list(structure_resnums)
    else:
        raise ValueError('Input must either be a list of residue numbers or a dictionary of residue numbers '
                         'and their frequency.')

    colors = sns.color_palette("hls", len(structure_resnums)).as_hex()

    to_show_chains = '( '
    for c in chain_id:
        to_show_chains += ':{} or'.format(c)
    to_show_chains = to_show_chains.strip(' or ')
    to_show_chains += ' )'

    for i, x in enumerate(structure_resnums):
        if isinstance(x, tuple):
            to_show_res = '( '
            for mut in x:
                to_show_res += '{} or '.format(mut)
            to_show_res = to_show_res.strip(' or ')
            to_show_res += ' )'
        else:
            to_show_res = x

        log.info('Selection: {} and not hydrogen and {}'.format(to_show_chains, to_show_res))

        if unique_colors:
            view.add_ball_and_stick(selection='{} and not hydrogen and {}'.format(to_show_chains, to_show_res),
                                    color=colors[i], opacity=opacity_dict[x], scale=scale_dict[x])
        else:
            view.add_ball_and_stick(selection='{} and not hydrogen and {}'.format(to_show_chains, to_show_res),
                                    color=color, opacity=opacity_dict[x], scale=scale_dict[x])


def add_features_to_nglview(view, structure_resnums, chain_id):
    """Add select features from the selected SeqProp object to an NGLWidget view object.

    Currently parsing for:
    * Single residue features (ie. metal binding sites)
    * Disulfide bonds

    Args:
        view (NGLWidget): NGLWidget view object
        seqprop (SeqProp): SeqProp object
        structprop (StructProp): StructProp object
        chain_id (str): ID of the structure's chain to get annotation from

    """
    # Parse and store chain seq if not already stored
    if not structprop.chains.has_id(chain_id):
        structprop.parse_structure()
        if not structprop.chains.has_id(chain_id):
            raise ValueError('Chain {} not present in structure {}'.format(chain_id, structprop.id))

    if not seqprop.features:
        log.warning('{}: no stored features'.format(seqprop.id))

    # Loop through any stored features
    for f in seqprop.features:

        # Display disulfide bonds
        if f.type.lower() == 'disulfide bond':
            # TODO: double check if .start or .start + 1
            disulfide = map_seqprop_resnums_to_structprop_resnums(resnums=[f.location.start + 1, f.location.end],
                                                                  seqprop=seqprop,
                                                                  structprop=structprop,
                                                                  chain_id=chain_id,
                                                                  use_representatives=False)
            to_view = [str(x)+'.CA' for x in list(disulfide.values())]
            view.add_distance(atom_pair=[to_view], color='black')
            log.info('Disulfide bridge at residues {} & {}'.format(f.location.start + 1, f.location.end))

        # Display DNA-binding regions
        if f.type.lower() == 'dna-binding region' or f.type.lower() == 'nucleotide phosphate-binding region':
            impres = self.map_seqprop_resnums_to_structprop_resnums(resnums=[f.location.start + 1,
                                                                             f.location.end],
                                                                    seqprop=seqprop,
                                                                    structprop=structprop,
                                                                    chain_id=chain_id,
                                                                    use_representatives=use_representatives)

            # TODO: need to check if f.location.start was mapped and if not, try incrementing. or  input the list
            # of resnums, not just the start and end
            if f.location.start + 1 in impres and f.location.end in impres:
                mapped_start = impres[f.location.start + 1]
                mapped_end = impres[f.location.end]
                view.add_ball_and_stick(selection=':{} and ( {}-{} )'.format(chain_id,
                                                                             mapped_start,
                                                                             mapped_end), color='black')
                log.info('{} at sequence region {}-{}, structure residues {}-{}'.format(f.type,
                                                                                        f.location.start,
                                                                                        f.location.end,
                                                                                        mapped_start,
                                                                                        mapped_end))

        # Display other single residues
        if f.location.end - 1 == f.location.start:
            if f.type.lower() == 'sequence variant' or f.type.lower() == 'mutagenesis site':
                continue
            impres = self.map_seqprop_resnums_to_structprop_resnums(resnums=f.location.end,
                                                                    seqprop=seqprop,
                                                                    structprop=structprop,
                                                                    chain_id=chain_id,
                                                                    use_representatives=use_representatives)
            if f.location.end in impres:
                impres_mapped = impres[f.location.end]
                view.add_ball_and_stick(selection=str(impres_mapped), color='black')
                view.add_label(selection=':{} and {}'.format(chain_id, impres_mapped), label_type='res', color='black')
                log.info('{} at sequence residue {}, structure residue {}'.format(f.type, f.location.end, impres_mapped))

def add_mutations_to_nglview(self, view, alignment_type='seqalign', alignment_ids=None,
                             seqprop=None, structprop=None, chain_id=None, use_representatives=False,
                             grouped=False, color='red', unique_colors=True,
                             opacity_range=(0.8,1), scale_range=(1,5)):
    """Add representations to an NGLWidget view object for residues that are mutated in the
        ``sequence_alignments`` attribute.

    Args:
        view (NGLWidget): NGLWidget view object
        alignment_type (str): Specified alignment type contained in the ``annotation`` field of an alignment object,
            ``seqalign`` or ``structalign`` are the current types.
        alignment_ids (str, list): Specified alignment ID or IDs to use
        seqprop (SeqProp): SeqProp object
        structprop (StructProp): StructProp object
        chain_id (str): ID of the structure's chain to get annotation from
        use_representatives (bool): If the representative sequence/structure/chain IDs should be used
        grouped (bool): If groups of mutations should be colored and sized together
        color (str): Color of the mutations (overridden if unique_colors=True)
        unique_colors (bool): If each mutation/mutation group should be colored uniquely
        opacity_range (tuple): Min/max opacity values (mutations that show up more will be opaque)
        scale_range (tuple): Min/max size values (mutations that show up more will be bigger)

    """
    if use_representatives:
        if seqprop and structprop and chain_id:
            raise ValueError('Overriding sequence, structure, and chain IDs with representatives. '
                             'Set use_representatives to False if custom IDs are to be used.')
    else:
        if not seqprop or not structprop or not chain_id:
            raise ValueError('Input sequence, structure, and chain to map between, or set use_representatives '
                             'to True.')

    if use_representatives:
        seqprop = self.representative_sequence
        structprop = self.representative_structure
        chain_id = self.representative_chain

    log.debug('Using sequence: {}, structure: {}, chain: {}'.format(seqprop.id, structprop.id, chain_id))

    # Get the summary of mutations
    single, fingerprint = self.sequence_mutation_summary(alignment_type=alignment_type, alignment_ids=alignment_ids)

    # Map residues from sequence to structure
    if not grouped:
        single_lens = {k: len(v) for k, v in single.items()}
        single_map_to_structure = {}

        for k, v in single_lens.items():
            resnum = int(k[1])
            resnum_to_structure = self.map_seqprop_resnums_to_structprop_resnums(resnums=resnum,
                                                                                 seqprop=seqprop,
                                                                                 structprop=structprop,
                                                                                 chain_id=chain_id,
                                                                                 use_representatives=use_representatives)
            if resnum not in resnum_to_structure:
                log.warning('{}: residue is not available in structure {}'.format(resnum, structprop.id))
                continue
            new_key = resnum_to_structure[resnum]
            single_map_to_structure[new_key] = v

        structprop.add_scaled_residues_highlight_to_nglview(view=view,
                                                            structure_resnums=single_map_to_structure,
                                                            chain=chain_id,
                                                            color=color,
                                                            unique_colors=unique_colors,
                                                            opacity_range=opacity_range,
                                                            scale_range=scale_range)
    else:
        log.warning('Viewing mutation groups is currently in beta -- groups may overwrite each other')
        fingerprint_lens = {k: len(v) for k, v in fingerprint.items()}
        fingerprint_map_to_structure = {}

        for k, v in fingerprint_lens.items():
            k_list = [int(x[1]) for x in k]
            resnums_to_structure = self.map_seqprop_resnums_to_structprop_resnums(resnums=k_list,
                                                                                  seqprop=seqprop,
                                                                                  structprop=structprop,
                                                                                  chain_id=chain_id,
                                                                                  use_representatives=use_representatives)
            new_key = tuple(y for y in resnums_to_structure.values())
            fingerprint_map_to_structure[new_key] = v

        structprop.add_scaled_residues_highlight_to_nglview(view=view,
                                                            structure_resnums=fingerprint_map_to_structure,
                                                            chain=chain_id,
                                                            color=color,
                                                            unique_colors=unique_colors,
                                                            opacity_range=opacity_range,
                                                            scale_range=scale_range)

def add_fingerprint_to_nglview(self, view, fingerprint,
                               seqprop=None, structprop=None, chain_id=None, use_representatives=False,
                               color='red', opacity_range=(0.8, 1), scale_range=(1, 5)):
    """Add representations to an NGLWidget view object for residues that are mutated in the
        ``sequence_alignments`` attribute.

    Args:
        view (NGLWidget): NGLWidget view object
        fingerprint (dict): Single mutation group from the ``sequence_mutation_summary`` function
        seqprop (SeqProp): SeqProp object
        structprop (StructProp): StructProp object
        chain_id (str): ID of the structure's chain to get annotation from
        use_representatives (bool): If the representative sequence/structure/chain IDs should be used
        color (str): Color of the mutations (overridden if unique_colors=True)
        opacity_range (tuple): Min/max opacity values (mutations that show up more will be opaque)
        scale_range (tuple): Min/max size values (mutations that show up more will be bigger)

    """
    if use_representatives:
        if seqprop and structprop and chain_id:
            raise ValueError('Overriding sequence, structure, and chain IDs with representatives. '
                             'Set use_representatives to False if custom IDs are to be used.')
    else:
        if not seqprop or not structprop or not chain_id:
            raise ValueError('Input sequence, structure, and chain to map between, or set use_representatives '
                             'to True.')

    if use_representatives:
        seqprop = self.representative_sequence
        structprop = self.representative_structure
        chain_id = self.representative_chain

    log.debug('Using sequence: {}, structure: {}, chain: {}'.format(seqprop.id, structprop.id, chain_id))

    fingerprint_lens = {k: len(v) for k, v in fingerprint.items()}
    fingerprint_map_to_structure = {}
    for k, v in fingerprint_lens.items():
        k_list = [int(x[1]) for x in k]
        resnums_to_structure = self.map_seqprop_resnums_to_structprop_resnums(resnums=k_list,
                                                                              seqprop=seqprop,
                                                                              structprop=structprop,
                                                                              chain_id=chain_id,
                                                                              use_representatives=use_representatives)
        new_key = tuple(y for y in resnums_to_structure.values())
        fingerprint_map_to_structure[new_key] = v

    structprop.add_scaled_residues_highlight_to_nglview(view=view,
                                                        structure_resnums=fingerprint_map_to_structure,
                                                        chain=chain_id,
                                                        color=color,
                                                        opacity_range=opacity_range,
                                                        scale_range=scale_range)