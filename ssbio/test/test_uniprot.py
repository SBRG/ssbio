import unittest
import six

import ssbio.databases.uniprot
from ssbio.databases.uniprot import UniProtProp

class TestUniProt(unittest.TestCase):
    """Unit tests for UniProt
    """

    def test_uniprotprop(self):
        up = UniProtProp(ident='P0ABP8', sequence_path='test_files/sequences/P0ABP8.fasta',
                         metadata_path='test_files/sequences/P0ABP8.txt')

        up_keys = ['bigg','description','ec_number','entry_version','gene_name','id','kegg',
                   'metadata_file','metadata_path','pdbs','pfam','refseq','reviewed','sequence_file','seq_len',
                   'sequence_path','seq_version','uniprot']

        for k in up_keys:
            # print(k)
            self.assertTrue(hasattr(up, k))

        self.assertEqual(up.bigg, None)
        self.assertEqual(up.description, ['Inosine phosphorylase', 'PNP', 'Purine nucleoside phosphorylase DeoD-type'])
        self.assertEqual(up.ec_number, ['2.4.2.1'])
        self.assertEqual(up.entry_version, '2016-11-02')
        self.assertEqual(up.gene_name, 'deoD')
        self.assertEqual(up.id, 'P0ABP8')
        six.assertCountEqual(self, up.kegg, ['ecj:JW4347', 'eco:b4384'])
        self.assertEqual(up.metadata_file, 'P0ABP8.txt')
        six.assertCountEqual(self, up.pdbs,
                        ['1OUM', '4TTA', '4TS3', '1OTY', '1K9S', '3OOE', '1OVG', '3ONV', '3UT6', '1OV6', '4TTJ', '1OTX',
                         '1OU4', '3OOH', '4TS9', '1A69', '3OPV', '1ECP', '4TTI'])
        self.assertEqual(up.pfam, ['PF01048'])
        six.assertCountEqual(self, up.refseq, ['NP_418801.1', 'NC_000913.3', 'WP_000224877.1', 'NZ_LN832404.1'])
        self.assertEqual(up.reviewed, True)
        self.assertEqual(up.sequence_file, 'P0ABP8.fasta')
        self.assertEqual(up.seq_len, 239)
        self.assertEqual(up.seq_version, '2007-01-23')
        self.assertEqual(up.uniprot, 'P0ABP8')

    def test_uniprot_valid_id(self):
        valid_ids = ['P20020',
                     'V9HWB9',
                     'A0A024QZK7',
                     'Q9Y633',
                     'O75414',
                     'B7Z3W8',
                     'V9HWJ2',
                     'F2Z2Y4',
                     'A0A024RDQ9',
                     'P48426']
        invalid_ids = ['1AAG', '1AN3', '1A8X', '8CEL', 'ASDASDASD##']

        for u in valid_ids:
            self.assertTrue(ssbio.databases.uniprot.is_valid_uniprot_id(u))

        for u in invalid_ids:
            self.assertFalse(ssbio.databases.uniprot.is_valid_uniprot_id(u))

    def test_uniprot_reviewed_checker(self):
        status = {'I1Z9G4': False,
                  'O14735': True}

        for f,v in status.items():
            self.assertEqual(ssbio.databases.uniprot.uniprot_reviewed_checker(f), v)

    def test_uniprot_reviewed_checker_batch(self):
        status = {'I1Z9G4': False,
                  'I3L3I9': False,
                  'J3KN10': False,
                  'J3KQV8': False,
                  'J3KRC4': False,
                  'O00408': True,
                  'O00757': True,
                  'O00764': True,
                  'O14494': True,
                  'O14495': True,
                  'O14735': True}

        test_status = ssbio.databases.uniprot.uniprot_reviewed_checker_batch(list(status.keys()))
        self.assertEqual(test_status, status)

    # def test_uniprot_metadata(self):
    #     test_out = {
    #         'P36959': {'u_description': ['Guanosine monophosphate reductase 1 {ECO:0000255|HAMAP-Rule:MF_03195}',
    #                                      "Guanosine 5'-monophosphate oxidoreductase 1 {ECO:0000255|HAMAP-Rule:MF_03195}",
    #                                      'GMPR 1 {ECO:0000255|HAMAP-Rule:MF_03195}',
    #                                      'GMP reductase 1 {ECO:0000255|HAMAP-Rule:MF_03195}'],
    #                    'u_ec_number': ['1.7.1.7'],
    #                    'u_entry_version': '2016-04-13',
    #                    'u_gene_name': 'GMPR',
    #                    'u_go': ['GO:0005829; C:cytosol; TAS:Reactome.',
    #                             'GO:1902560; C:GMP reductase complex; IEA:UniProtKB-EC.',
    #                             'GO:0003920; F:GMP reductase activity; TAS:ProtInc.',
    #                             'GO:0046872; F:metal ion binding; IEA:UniProtKB-KW.',
    #                             'GO:0055086; P:nucleobase-containing small molecule metabolic process; TAS:Reactome.',
    #                             'GO:0006144; P:purine nucleobase metabolic process; TAS:Reactome.',
    #                             'GO:0006163; P:purine nucleotide metabolic process; IEA:UniProtKB-HAMAP.',
    #                             'GO:0043101; P:purine-containing compound salvage; TAS:Reactome.',
    #                             'GO:0009409; P:response to cold; TAS:ProtInc.',
    #                             'GO:0044281; P:small molecule metabolic process; TAS:Reactome.'],
    #                    'u_kegg_id': ['hsa:2766'],
    #                    'u_pfam': ['PF00478'],
    #                    'u_refseq': ['NP_006868.3', 'NM_006877.3'],
    #                    'u_reviewed': True,
    #                    'u_seq': 'MPRIDADLKLDFKDVLLRPKRSSLKSRAEVDLERTFTFRNSKQTYSGIPIIVANMDTVGTFEMAAVMSQHSMFTAIHKHYSLDDWKLFATNHPECLQNVAVSSGSGQNDLEKMTSILEAVPQVKFICLDVANGYSEHFVEFVKLVRAKFPEHTIMAGNVVTGEMVEELILSGADIIKVGVGPGSVCTTRTKTGVGYPQLSAVIECADSAHGLKGHIISDGGCTCPGDVAKAFGAGADFVMLGGMFSGHTECAGEVFERNGRKLKLFYGMSSDTAMNKHAGGVAEYRASEGKTVEVPYKGDVENTILDILGGLRSTCTYVGAAKLKELSRRATFIRVTQQHNTVFS',
    #                    'u_seq_len': 345,
    #                    'u_seq_version': '1994-06-01',
    #                    'u_uniprot_acc': 'P36959'}}
    #
    #     self.assertEqual(ssbio.databases.uniprot.uniprot_metadata('P36959'), test_out)