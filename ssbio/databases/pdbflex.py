import requests
import ssbio.utils
import os.path as op


# #### PDB stats
# Request flexibility data about one particular PDB.
#
# http://pdbflex.org/php/api/PDBStats.php?pdbID=1a50&chainID=A
#
#     pdbID   of structure you are interested in
#     chainID   of chain you are interested in
#
#     [{"pdbID":"1a50",
#     "chainID":"A",
#     "parentClusterID":"4hn4A",
#     "avgRMSD":"0.538",
#     "maxRMSD":"2.616",
#     "flexibilityLabel":"Low",
#     "otherClusterMembers":["4hn4A","4hpjA","4hpxA","4kkxA",...],
#     "PDBFlexLink":"http:\/\/pdbflex.org\/cluster.html#!\/4hn4A\/20987\/1a50A"}]
#
# Note: you can omit the chainID and PDBFlex will return information for all chains.
#
# #### RMSD profile
# Request RMSD array used for local flexibility plots
#
# http://pdbflex.org/php/api/rmsdProfile.php?pdbID=1a50&chainID=A
#
# pdbID   PDB ID of structure you are interested in
# chainID   Chain ID of chain you are interested in
#
#     {"queryPDB":"1a50A",
#     "clusterName":"4hn4A",
#     "profile":"[0.616,0.624,0.624,0.624,0.624,0.624,0.029,0.013,0.016,0.023,0.025,0.028,0.030,0.034,0.035,0.035,0.035,0.035,0.036,0.033,0.027,0.023,0.017...]"}
#
# #### PDB representatives
# Request representatives for a PDB's own cluster. Returns a list of chains that represent the most distinct structures in the cluster.
#
# http://pdbflex.org/php/api/representatives.php?pdbID=1a50&chainID=A
#
#     pdbID   PDB ID of structure you are interested in
#     chainID   Chain ID of chain you are interested in
#
#     ["2trsA","3pr2A","1kfjA"]


def get_pdbflex_info(pdb_id, chain_id, outdir, force_rerun=False):
    outfile = '{}{}_pdbflex_stats.json'.format(pdb_id, chain_id)

    pdbflex_link = 'http://pdbflex.org/php/api/PDBStats.php?pdbID={}&chainID={}'.format(pdb_id,
                                                                                        chain_id)
    infolist = ssbio.utils.request_json(link=pdbflex_link, outfile=outfile, outdir=outdir, force_rerun_flag=force_rerun)

    # TODO: will running with chain ID always return a single item list?
    assert len(infolist) == 1

    newdict = {}
    for k, v in infolist[0].items():
        if k == 'avgRMSD' and v:
            newdict[k] = float(v)
        elif k == 'maxRMSD' and v:
            newdict[k] = float(v)
        else:
            newdict[k] = v

    return newdict


def get_pdbflex_rmsd_profile(pdb_id, chain_id, outdir, force_rerun=False):
    outfile = '{}{}_pdbflex_rmsdprofile.json'.format(pdb_id, chain_id)

    pdbflex_link = 'http://pdbflex.org/php/api/rmsdProfile.php?pdbID={}&chainID={}'.format(pdb_id,
                                                                                           chain_id)
    infodict = ssbio.utils.request_json(link=pdbflex_link, outfile=outfile, outdir=outdir, force_rerun_flag=force_rerun)
    infodict['profile'] = [float(x) for x in infodict['profile'].strip('[]').split(',')]
    return infodict


def get_pdbflex_representatives(pdb_id, chain_id, outdir, force_rerun=False):
    outfile = '{}{}_pdbflex_representatives.json'.format(pdb_id, chain_id)

    pdbflex_link = 'http://pdbflex.org/php/api/representatives.php?pdbID={}&chainID={}'.format(pdb_id,
                                                                                               chain_id)
    infolist = ssbio.utils.request_json(link=pdbflex_link, outfile=outfile, outdir=outdir, force_rerun_flag=force_rerun)
    #     infolist = [str(x) for x in infolist.strip('[]').split(',')]
    return infolist