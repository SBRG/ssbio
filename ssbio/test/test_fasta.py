import os
import os.path as op
import tempfile
import unittest
import filecmp
import ssbio.sequence.utils.fasta as fasta


class TestFasta(unittest.TestCase):
    """Unit tests for sequence/fasta tools."""

    def test_load_fasta_file(self):
        sequences = fasta.load_fasta_file('test_sequences/X5D299-1.faa')
        self.assertEqual(len(sequences), 1)
        self.assertEqual(sequences[0].id, 'X5D299-1')
        self.assertRaises(IOError, fasta.load_fasta_file, 'NOTAFILE.faa')

    def test_load_seq_string_as_seqrecord(self):
        test_string = 'ASQAGIPSGVYNVIPCSRKNAKEVGEAICTDPLVSKISF'
        ident = 'TEST'
        sequence = fasta.load_seq_string_as_seqrecord(test_string, ident)
        self.assertEqual(len(sequence), len(test_string))
        self.assertEqual(sequence.id, ident)

        self.assertRaises(ValueError, fasta.load_seq_string_as_seqrecord, 'N-0t,V4L1DAA', ident)

    def test_write_fasta_file(self):
        test_string = 'ASQAGIPSGVYNVIPCSRKNAKEVGEAICTDPLVSKISF'
        ident = 'TEST'
        sequence = fasta.load_seq_string_as_seqrecord(test_string, ident, desc='tester')
        custom_ext = '.fasta'
        custom_outname = 'TESTER'
        overwrite = True

        custom_dir = op.join('test_files', 'out')
        gold_standard = op.join('test_files', 'sequences', 'TESTER.faa')

        written_file_cusdir = fasta.write_fasta_file(sequence, custom_outname,
                                                     outdir=custom_dir, force_rerun=overwrite)
        self.assertEqual(written_file_cusdir, os.path.join(custom_dir, custom_outname + '.faa'))
        # Test is file contents are equal with filecmp
        self.assertTrue(filecmp.cmp(gold_standard, written_file_cusdir))

        written_file_cusdirid = fasta.write_fasta_file(sequence, custom_outname,
                                                       outdir=custom_dir, force_rerun=overwrite)
        self.assertEqual(written_file_cusdirid, os.path.join(custom_dir, custom_outname + '.faa'))
        self.assertTrue(filecmp.cmp(gold_standard, written_file_cusdirid))

        written_file_cusdiridext = fasta.write_fasta_file(sequence, custom_outname,
                                                          outdir=custom_dir, outext=custom_ext, force_rerun=overwrite)
        self.assertEqual(written_file_cusdiridext, os.path.join(custom_dir, custom_outname + custom_ext))
        self.assertTrue(filecmp.cmp(gold_standard, written_file_cusdiridext))

    def test_write_fasta_file_from_dict(self):
        tester = {'TEST': 'ASQAGIPSGVYNVIPCSRKNAKEVGEAICTDPLVSKISF'}
        custom_ext = '.fasta'
        custom_outname = 'TESTER'
        overwrite = True

        with tempfile.TemporaryDirectory() as custom_dir:
            written_file_cusdir = fasta.write_fasta_file_from_dict(tester, custom_outname,
                                                                   outdir=custom_dir, force_rerun=overwrite)
            self.assertEqual(written_file_cusdir, os.path.join(custom_dir, custom_outname + '.faa'))

            written_file_cusdirid = fasta.write_fasta_file_from_dict(tester, custom_outname,
                                                                     outdir=custom_dir, force_rerun=overwrite)
            self.assertEqual(written_file_cusdirid, os.path.join(custom_dir, custom_outname + '.faa'))

            written_file_cusdiridext = fasta.write_fasta_file_from_dict(tester, custom_outname,
                                                                        outdir=custom_dir, outext=custom_ext, force_rerun=overwrite)
            self.assertEqual(written_file_cusdiridext, os.path.join(custom_dir, custom_outname + custom_ext))

if __name__ == "__main__":
    unittest.main()
