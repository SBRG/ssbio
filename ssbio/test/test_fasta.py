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
        tester = {'TEST': 'ASQAGIPSGVYNVIPCSRKNAKEVGEAICTDPLVSKISF'}
        custom_ext = '.fasta'
        custom_outname = 'TESTER'
        overwrite = True

        with tempfile.TemporaryDirectory() as custom_dir:
            written_file_cusdir = fasta.write_fasta_file(tester, custom_outname,
                                                         outdir=custom_dir, force_rerun=overwrite)
            self.assertEqual(written_file_cusdir, os.path.join(custom_dir, custom_outname + '.faa'))

            written_file_cusdirid = fasta.write_fasta_file(tester, custom_outname,
                                                           outdir=custom_dir, force_rerun=overwrite)
            self.assertEqual(written_file_cusdirid, os.path.join(custom_dir, custom_outname + '.faa'))

            written_file_cusdiridext = fasta.write_fasta_file(tester, custom_outname,
                                                              outdir=custom_dir, outext=custom_ext, force_rerun=overwrite)
            self.assertEqual(written_file_cusdiridext, os.path.join(custom_dir, custom_outname + custom_ext))

if __name__ == "__main__":
    unittest.main()
