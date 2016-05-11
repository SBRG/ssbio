import unittest
import os
import tempfile
import ssbio.sequence.fasta as fasta


class TestFasta(unittest.TestCase):
    """Unit tests for sequence/fasta tools."""

    def test_load_fasta_file(self):
        sequences = fasta.load_fasta_file('test_sequences/X5D299-1.faa')
        self.assertEqual(len(sequences), 1)
        self.assertEqual(sequences[0].id, 'X5D299-1')
        self.assertRaises(IOError, fasta.load_fasta_file, 'NOTAFILE.faa')

    def test_load_seq_string(self):
        test_string = 'ASQAGIPSGVYNVIPCSRKNAKEVGEAICTDPLVSKISF'
        ident = 'TEST'
        sequence = fasta.load_seq_string(test_string, ident)
        self.assertEqual(len(sequence), len(test_string))
        self.assertEqual(sequence.id, ident)
        self.assertRaises(ValueError, fasta.load_seq_string, 'N-0t,V4L1DAA', ident)

    def test_write_fasta_file(self):
        sequence = 'ASQAGIPSGVYNVIPCSRKNAKEVGEAICTDPLVSKISF'
        ident = 'TEST'
        custom_ext = 'fasta'
        custom_id = 'TEST2'
        overwrite = True

        with tempfile.TemporaryDirectory() as custom_dir:
            written_file_cusdir = fasta.write_fasta_file(
                sequence, ident, outpath=custom_dir, overwrite=overwrite)
            self.assertEqual(written_file_cusdir, os.path.join(
                custom_dir, ident + '.faa'))

            written_file_cusdirid = fasta.write_fasta_file(
                sequence, custom_id, outpath=custom_dir, overwrite=overwrite)
            self.assertEqual(written_file_cusdirid, os.path.join(
                custom_dir, custom_id + '.faa'))

            written_file_cusdiridext = fasta.write_fasta_file(
                sequence, custom_id, outpath=custom_dir, extension=custom_ext, overwrite=overwrite)
            self.assertEqual(written_file_cusdiridext, os.path.join(
                custom_dir, custom_id + '.' + custom_ext))

            # TODO: test if file was overwritten? how to do that?
            # written_file_cusdirnew = s.write_fasta_file(sequence, outdir=custom_dir, overwrite=overwrite)
            # self.assertEqual(written_file_cusdirnew, os.path.join(custom_dir, ident + '.faa'))

if __name__ == "__main__":
    unittest.main()
