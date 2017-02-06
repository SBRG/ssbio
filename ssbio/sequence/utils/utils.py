from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein


def cast_to_str(obj):
    """Return a string representation of a Seq or SeqRecord

    Args:
        obj:

    Returns:

    """
    if isinstance(obj, str):
        return obj
    if isinstance(obj, Seq):
        return str(obj)
    if isinstance(obj, SeqRecord):
        return str(obj.seq)
    else:
        raise ValueError('Must provide a string, Seq, or SeqRecord object.')


def cast_to_seq_record(obj, alphabet=generic_protein, id="<unknown id>", name="<unknown name>",
                       description="<unknown description>", dbxrefs=None,
                       features=None, annotations=None,
                       letter_annotations=None):
    if isinstance(obj, SeqRecord):
        return obj
    if isinstance(obj, Seq):
        return SeqRecord(obj, id, name, description, dbxrefs, features, annotations, letter_annotations)
    if isinstance(obj, str):
        obj = obj.upper()
        return SeqRecord(Seq(obj, alphabet), id, name, description, dbxrefs, features, annotations, letter_annotations)
    else:
        raise ValueError('Must provide a string, Seq, or SeqRecord object.')