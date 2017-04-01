import unittest

from ssbio.protein.sequence.properties.scratch import SCRATCH


class TestSCRATCH(unittest.TestCase):
    """Unit tests for SCRATCH
    """

    @classmethod
    def setUpClass(self):
        try:
            self.my_scratch = SCRATCH(project_name='test', seq_file='test_files/scratch/test.fasta')
            self.my_scratch.run_scratch(
                path_to_scratch='/home/nathan/software/SCRATCH-1D_1.1/bin/run_SCRATCH-1D_predictors.sh',
                outname='test', outdir='test_files/scratch/', num_cores=4,
                force_rerun=False)
        except OSError:
            pass

    def test_sspro_summary(self):
        results = {'pdb,4i8h,A,#1': {'H': 0.1031390134529148,
                                     'E': 0.34080717488789236,
                                     'C': 0.5560538116591929},
                   'pdb,3w5h,A,#1': {'H': 0.21691176470588236,
                                     'E': 0.35294117647058826,
                                     'C': 0.43014705882352944},
                   'pdb,4i8g,A,#1': {'H': 0.1031390134529148,
                                     'E': 0.34080717488789236,
                                     'C': 0.5560538116591929},
                   'pdb,4eic,A,#1': {'H': 0.5591397849462365,
                                     'E': 0.021505376344086023,
                                     'C': 0.41935483870967744}}

        sspro_dict = self.my_scratch.sspro_summary()

        for k,v in results.items():
            self.assertAlmostEqual(results[k], sspro_dict[k])

    def test_sspro8_summary(self):
        results = {'pdb,4i8g,A,#1': {'H': 0.07174887892376682, 'C': 0.2556053811659193, 'I': 0.0,
                                     'E': 0.32286995515695066, 'G': 0.03139013452914798, 'T': 0.14798206278026907,
                                     'B': 0.017937219730941704, 'S': 0.15246636771300448},
                   'pdb,3w5h,A,#1': {'H': 0.17279411764705882, 'C': 0.21323529411764705, 'I': 0.022058823529411766,
                                     'E': 0.34558823529411764, 'G': 0.04411764705882353, 'T': 0.11397058823529412,
                                     'B': 0.007352941176470588, 'S': 0.08088235294117647},
                   'pdb,4eic,A,#1': {'H': 0.5268817204301075, 'C': 0.13978494623655913, 'I': 0.0, 'E': 0.0,
                                     'G': 0.03225806451612903, 'T': 0.17204301075268819, 'B': 0.021505376344086023,
                                     'S': 0.10752688172043011},
                   'pdb,4i8h,A,#1': {'H': 0.07174887892376682, 'C': 0.2556053811659193, 'I': 0.0,
                                     'E': 0.32286995515695066,
                                     'G': 0.03139013452914798, 'T': 0.14798206278026907, 'B': 0.017937219730941704,
                                     'S': 0.15246636771300448}}

        sspro_dict = self.my_scratch.sspro8_summary()

        for k,v in results.items():
            self.assertAlmostEqual(results[k], sspro_dict[k])

    def test_accpro_summary(self):
        results = {'pdb,4eic,A,#1': {'exposed': 0.6344086021505376, 'buried': 0.3655913978494624},
                   'pdb,4i8h,A,#1': {'exposed': 0.4439461883408072, 'buried': 0.5560538116591929},
                   'pdb,3w5h,A,#1': {'exposed': 0.5036764705882353, 'buried': 0.4963235294117647},
                   'pdb,4i8g,A,#1': {'exposed': 0.4439461883408072, 'buried': 0.5560538116591929}}

        accpro_dict = self.my_scratch.accpro_summary()

        for k,v in results.items():
            self.assertAlmostEqual(results[k], accpro_dict[k])

    def test_accpro20_summary(self):
        cutoff = 25
        results = {'pdb,4eic,A,#1': {'buried': 0.41935483870967744, 'exposed': 0.5806451612903226},
                   'pdb,4i8g,A,#1': {'buried': 0.6143497757847534, 'exposed': 0.38565022421524664},
                   'pdb,4i8h,A,#1': {'buried': 0.6143497757847534, 'exposed': 0.38565022421524664},
                   'pdb,3w5h,A,#1': {'buried': 0.5404411764705882, 'exposed': 0.45955882352941174}}

        accpro_dict = self.my_scratch.accpro20_summary(cutoff=cutoff)

        for k,v in results.items():
            self.assertAlmostEqual(results[k], accpro_dict[k])

    @classmethod
    def tearDownClass(self):
        pass