import pytest
import os.path as op
import ssbio.databases.pdb as pdb
from six.moves.urllib_error import URLError


def test_download_mmcif_header(pdb_ids_working, pdb_ids_obsolete, pdb_ids_false, test_files_tempdir):
    for wp in pdb_ids_working:
        pdb.download_mmcif_header(pdb_id=wp, outdir=test_files_tempdir, force_rerun=True)
        outfile = '{}.header.cif'.format(wp)
        assert op.isfile(op.join(test_files_tempdir, outfile))
    # for obp in pdb_ids_obsolete:  # TODO: you can still download obsolete headers actually
    #     with pytest.raises(URLError):
    #         pdb.download_mmcif_header(pdb_id=obp, outdir=test_files_tempdir, force_rerun=True)
    for fp in pdb_ids_false:
        with pytest.raises(URLError):
            pdb.download_mmcif_header(pdb_id=fp, outdir=test_files_tempdir, force_rerun=True)

def test_download_sifts_xml(pdb_ids_working, pdb_ids_obsolete, pdb_ids_false, test_files_tempdir):
    for wp in pdb_ids_working:
        pdb.download_sifts_xml(pdb_id=wp, outdir=test_files_tempdir)
        outfile = '{}.sifts.xml'.format(wp)
        assert op.isfile(op.join(test_files_tempdir, outfile))
    # for obp in pdb_ids_obsolete:  # TODO: you can still download obsolete sifts actually
    #     with pytest.raises(URLError):
    #         pdb.download_sifts_xml(pdb_id=obp, outdir=test_files_tempdir, force_rerun=True)
    for fp in pdb_ids_false:
        with pytest.raises(URLError):
            pdb.download_sifts_xml(pdb_id=fp, outdir=test_files_tempdir, force_rerun=True)

def test_map_uniprot_resnum_to_pdb(pdb_ids_working, pdb_ids_obsolete, pdb_ids_false, test_files_tempdir):
    pass

def test_best_structures(pdb_ids_working, pdb_ids_obsolete, pdb_ids_false, test_files_tempdir):
    pass

def test_blast_pdb(pdb_ids_working, pdb_ids_obsolete, pdb_ids_false, test_files_tempdir):
    pass

def test_blast_pdb_df(pdb_ids_working, pdb_ids_obsolete, pdb_ids_false, test_files_tempdir):
    pass

def test_get_resolution(pdb_ids_working, pdb_ids_obsolete, pdb_ids_false, test_files_tempdir):
    pass

def test_get_release_date(pdb_ids_working, pdb_ids_obsolete, pdb_ids_false, test_files_tempdir):
    pass

def test_get_num_bioassemblies(pdb_ids_working, pdb_ids_obsolete, pdb_ids_false, test_files_tempdir):
    pass

def test_get_bioassembly_info(pdb_ids_working, pdb_ids_obsolete, pdb_ids_false, test_files_tempdir):
    pass

def test_download_biomol(pdb_ids_working, pdb_ids_obsolete, pdb_ids_false, test_files_tempdir):
    for wp in pdb_ids_working:
        pdb.download_biomol(pdb_id=wp, outdir=test_files_tempdir)
        outfile = '{}.sifts.xml'.format(wp)
        assert op.isfile(op.join(test_files_tempdir, outfile))
    # for obp in pdb_ids_obsolete:  # TODO: you can still download obsolete sifts actually
    #     with pytest.raises(URLError):
    #         pdb.download_biomol(pdb_id=obp, outdir=test_files_tempdir, force_rerun=True)
    for fp in pdb_ids_false:
        with pytest.raises(URLError):
            pdb.download_biomol(pdb_id=fp, outdir=test_files_tempdir, force_rerun=True)