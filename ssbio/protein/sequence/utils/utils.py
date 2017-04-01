from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


def cast_to_str(obj):
    """Return a string representation of a Seq or SeqRecord.

    Args:
        obj (str, Seq, SeqRecord): Biopython Seq or SeqRecord

    Returns:
        str: String representation of the sequence

    """

    if isinstance(obj, str):
        return obj
    if isinstance(obj, Seq):
        return str(obj)
    if isinstance(obj, SeqRecord):
        return str(obj.seq)
    else:
        raise ValueError('Must provide a string, Seq, or SeqRecord object.')


def cast_to_seq_record(obj, alphabet=IUPAC.extended_protein, id="<unknown id>", name="<unknown name>",
                       description="<unknown description>", dbxrefs=None,
                       features=None, annotations=None,
                       letter_annotations=None):
    """Return a SeqRecord representation of a string or Seq object.

    Args:
        obj (str, Seq, SeqRecord): Sequence string or Biopython Seq object
        alphabet: See Biopython SeqRecord docs
        id: See Biopython SeqRecord docs
        name: See Biopython SeqRecord docs
        description: See Biopython SeqRecord docs
        dbxrefs: See Biopython SeqRecord docs
        features: See Biopython SeqRecord docs
        annotations: See Biopython SeqRecord docs
        letter_annotations: See Biopython SeqRecord docs

    Returns:
        SeqRecord: SeqRecord representation of the sequence

    """

    if isinstance(obj, SeqRecord):
        return obj
    if isinstance(obj, Seq):
        return SeqRecord(obj, id, name, description, dbxrefs, features, annotations, letter_annotations)
    if isinstance(obj, str):
        obj = obj.upper()
        return SeqRecord(Seq(obj, alphabet), id, name, description, dbxrefs, features, annotations, letter_annotations)
    else:
        raise ValueError('Must provide a string, Seq, or SeqRecord object.')