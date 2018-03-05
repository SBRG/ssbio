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

    def get_subsequence_from_property(self, property_key, property_value, condition):
        """Get a subsequence as a new SeqProp object given a certain property you want to find in
        this chain's letter_annotation

        See documentation for :func:`ssbio.protein.sequence.seqprop.SeqProp.get_subsequence_from_property`

        Args:
            property_key (str): Property key in the ``letter_annotations`` attribute that you want to filter using
            property_value (str): Property value that you want to filter by
            condition (str): ``<``, ``=``, ``>``, ``>=``, or ``<=`` to filter the values by

        Returns:
            SeqProp: New SeqProp object that you can run computations on or just extract its properties

        """
        if not self.seq_record:
            raise ValueError('No chain sequence stored')

        if property_key not in self.seq_record.letter_annotations:
            raise KeyError('{}: {} not contained in the letter annotations'.format(self.seq_record.id, property_key))

        subfeat_indices = list(locate(self.seq_record.letter_annotations[property_key],
                                      lambda x: ssbio.utils.check_condition(x, condition, property_value)))

        biop_compound_list = []
        for idx in subfeat_indices:
            feat = FeatureLocation(idx, idx + 1)
            biop_compound_list.append(feat)

        sub_feature_location = CompoundLocation(biop_compound_list)
        sub_feature = sub_feature_location.extract(self.seq_record)

        new_sp = SeqProp(id='{}_{}_{}_{}_extracted'.format(self.id, property_key, condition, property_value),
                         seq=sub_feature)
        new_sp.letter_annotations = sub_feature.letter_annotations

        return new_sp