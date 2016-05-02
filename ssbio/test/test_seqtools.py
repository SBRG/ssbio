import unittest
import os
import tempfile
from ssbio.itasser import seqtools

s = seqtools.SeqTools()


class Test(unittest.TestCase):
    """Unit tests for seqtools."""

    def test_load_fasta_file(self):
        sequences = s.load_fasta_file('sequences/X5D299-1.faa')
        self.assertEqual(len(sequences), 1)
        self.assertEqual(sequences[0].id, 'X5D299-1')
        self.assertRaises(IOError, s.load_fasta_file, 'NOTAFILE.faa')

    def test_load_seq_string(self):
        test_string = 'ASQAGIPSGVYNVIPCSRKNAKEVGEAICTDPLVSKISF'
        ident = 'TEST'
        sequence = s.load_seq_string(test_string, ident)
        self.assertEqual(len(sequence), len(test_string))
        self.assertEqual(sequence.id, ident)
        self.assertRaises(ValueError, s.load_seq_string, 'N-0t,V4L1DAA', ident)

    def test_write_fasta_file(self):
        test_string = 'ASQAGIPSGVYNVIPCSRKNAKEVGEAICTDPLVSKISF'
        ident = 'TEST'
        custom_ext = 'fasta'
        custom_id = 'TEST2'
        overwrite = True

        sequence = s.load_seq_string(test_string, ident)

        with tempfile.TemporaryDirectory() as custom_dir:
            written_file_cusdir = s.write_fasta_file(
                sequence, outdir=custom_dir)
            self.assertEqual(written_file_cusdir, os.path.join(
                custom_dir, ident + '.faa'))

            written_file_cusdirid = s.write_fasta_file(
                sequence, outdir=custom_dir, custom_id=custom_id)
            self.assertEqual(written_file_cusdirid, os.path.join(
                custom_dir, custom_id + '.faa'))

            written_file_cusdiridext = s.write_fasta_file(
                sequence, outdir=custom_dir, custom_id=custom_id, extension=custom_ext)
            self.assertEqual(written_file_cusdiridext, os.path.join(
                custom_dir, custom_id + '.' + custom_ext))

            # TODO: test if file was overwritten? how to do that?
            # written_file_cusdirnew = s.write_fasta_file(sequence, outdir=custom_dir, overwrite=overwrite)
            # self.assertEqual(written_file_cusdirnew, os.path.join(custom_dir, ident + '.faa'))

if __name__ == "__main__":
    unittest.main()
