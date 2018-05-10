from ssbio.core.object import Object
import logging
log = logging.getLogger(__name__)
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from more_itertools import locate
import ssbio.utils
from ssbio.protein.sequence.seqprop import SeqProp

class ChainProp(Object):
    """Class for protein structural properties of a specific chain"""

    def __init__(self, ident, pdb_parent, seq_record=None, description=None):
        Object.__init__(self, id=ident, description=description)

        self.pdb_parent = pdb_parent
        self.seq_record = seq_record

        if not self.description:
            self.description = 'Chain {} from PDB parent {}'.format(self.id, self.pdb_parent)

    def reset_seq_record(self):
        self.seq_record = None

    def get_subsequence(self, resnums, new_id=None, copy_letter_annotations=True):
        """Get a subsequence as a new SeqProp object given a list of residue numbers"""
        # XTODO: documentation

        if not self.seq_record:
            raise ValueError('No chain sequence stored')

        biop_compound_list = []
        for resnum in resnums:
            feat = FeatureLocation(resnum - 1, resnum)
            biop_compound_list.append(feat)

        if len(biop_compound_list) == 0:
            log.info('Zero length subsequences')
            return
        elif len(biop_compound_list) == 1:
            log.debug('Subsequence only one residue long')
            sub_feature_location = biop_compound_list[0]
        else:
            sub_feature_location = CompoundLocation(biop_compound_list)

        sub_feature = sub_feature_location.extract(self.seq_record)

        if not new_id:
            new_id = '{}_subseq'.format(self.id)

        new_sp = SeqProp(id=new_id, seq=sub_feature)
        if copy_letter_annotations:
            new_sp.letter_annotations = sub_feature.letter_annotations
        return new_sp

    def get_subsequence_from_property(self, property_key, property_value, condition, return_resnums=False,
                                      copy_letter_annotations=True):
        """Get a subsequence as a new SeqProp object given a certain property you want to find in
        this chain's letter_annotation

        See documentation for :func:`ssbio.protein.sequence.seqprop.SeqProp.get_subsequence_from_property`

        Args:
            property_key (str): Property key in the ``letter_annotations`` attribute that you want to filter using
            property_value (str): Property value that you want to filter by
            condition (str): ``<``, ``=``, ``>``, ``>=``, or ``<=`` to filter the values by
            return_resnums (bool): If resnums should be returned as well

        Returns:
            SeqProp: New SeqProp object that you can run computations on or just extract its properties

        """
        if not self.seq_record:
            raise ValueError('No chain sequence stored')

        if property_key not in self.seq_record.letter_annotations:
            log.error(KeyError('{}: {} not contained in the letter annotations'.format(self.seq_record.id, property_key)))
            return

        if condition == 'in':
            subfeat_indices = list(locate(self.seq_record.letter_annotations[property_key],
                                          lambda x: x in property_value))
        else:
            subfeat_indices = list(locate(self.seq_record.letter_annotations[property_key],
                                          lambda x: ssbio.utils.check_condition(x, condition, property_value)))
        subfeat_resnums = [x + 1 for x in subfeat_indices]

        new_sp = self.get_subsequence(resnums=subfeat_resnums, new_id='{}_{}_{}_{}_extracted'.format(self.pdb_parent,
                                                                                                     self.id,
                                                                                                     property_key,
                                                                                                     condition,
                                                                                                     property_value),
                                      copy_letter_annotations=copy_letter_annotations)

        if return_resnums:
            return new_sp, subfeat_resnums
        else:
            return new_sp