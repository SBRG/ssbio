import pandas as pd

def parse_init_dat(infile):
    """Parse the main init.dat file which contains the modeling results

    The first line of the file init.dat contains stuff like:
    "120 easy  40   8"
    # TODO: what is this?

    Args:
        infile:

    Returns:

    """
    return

def parse_bsites_inf(infile):
    """Parse the Bsites.inf output file of COACH and returns a dataframe

    Bsites.inf contains the summary of COACH clustering results after all
        other prediction algorithms have finished
    For each site (cluster), there are three lines:
    Line 1: site number, c-score of coach prediction, cluster size
    Line 2: algorithm, PDB ID, ligand ID, center of binding site (cartesian coordinates),
        c-score of the algorithm’kegg prediction, binding residues from single template
    Line 3: Statistics of ligands in the cluster

    Args:
        infile (str): Path to Bsites.inf

    Returns:
        bsites_inf_df (Pandas DataFrame): Organized dataframe of results, columns defined below
            site_num: cluster which is the consensus binding site
            c_score: confidence score of the cluster prediction
            cluster_size: number of predictions within this cluster
            algorithm: main? algorithm used to make the prediction
            pdb_template_id: PDB ID of the template used to make the prediction
            pdb_template_chain: chain of the PDB which has the ligand
            pdb_ligand: predicted ligand to bind
            binding_location_coords: centroid of the predicted ligand position in the homology model
            c_score_method: confidence score for the main algorithm
            binding_residues: predicted residues to bind the ligand
            ligand_cluster_counts: number of predictions per ligand

    """

    bsites_inf_df = pd.DataFrame()

    with open(infile) as pp:
        lines = list(filter(None, (line.rstrip() for line in pp)))

        for i in range(len(lines) // 3):
            bsites_site_dict = {}

            line1 = lines[i * 3].split('\t')
            line2 = lines[i * 3 + 1].split('\t')
            line3 = lines[i * 3 + 2]

            bsites_site_dict['site_num'] = line1[0]
            bsites_site_dict['c_score'] = float(line1[1])
            bsites_site_dict['cluster_size'] = line1[2]

            bsites_site_dict['algorithm'] = line2[0]
            bsites_site_dict['pdb_template_id'] = line2[1][:4]
            bsites_site_dict['pdb_template_chain'] = line2[1][4]
            bsites_site_dict['pdb_ligand'] = line2[2]
            bsites_site_dict['binding_location_coords'] = tuple(float(x) for x in line2[3].split())

            # TODO: what'kegg the difference between this c-score and the cluster'kegg c-score?
            # how is the cluster'kegg c-score computed? it'kegg not the average c-score of all methods
            # also why are some COFACTOR c-scores >1?
            # 160411 - seems like the COFACTOR "BS-score" is being reported here, not it'kegg c-score...
            tmp_split = line2[4].split(' :')
            bsites_site_dict['c_score_method'] = tmp_split[0]
            bsites_site_dict['binding_residues'] = tmp_split[1]

            bsites_site_dict['ligand_cluster_counts'] = line3

            bsites_inf_df = bsites_inf_df.append(bsites_site_dict, ignore_index=True)

    if len(bsites_inf_df) == 0:
        return bsites_inf_df

    bsites_inf_df = bsites_inf_df[['site_num', 'c_score', 'cluster_size', 'algorithm',
                                   'pdb_template_id', 'pdb_template_chain', 'pdb_ligand',
                                   'binding_location_coords', 'c_score_method', 'binding_residues',
                                   'ligand_cluster_counts']]
    return bsites_inf_df  # .set_index('site_num')


def parse_ec(infile):
    """Parse the EC.dat output file of COACH and returns a dataframe

    EC.dat contains the predicted EC number and active residues.
        The columns are: PDB_ID, TM-score, RMSD, Sequence identity,
        Coverage, Confidence score, EC number, and Active site residues

    Args:
        infile (str): Path to EC.dat

    Returns:
        ec_df (Pandas DataFrame): Organized dataframe of results, columns defined below
            pdb_template_id: PDB ID of the template used to make the prediction
            pdb_template_chain: chain of the PDB which has the ligand
            tm_score: TM-score of the template to the model (similarity score)
            rmsd: RMSD of the template to the model (also a measure of similarity)
            seq_ident: percent sequence identity
            seq_coverage: percent sequence coverage
            c_score: confidence score of the EC prediction
            ec_number: predicted EC number
            binding_residues: predicted residues to bind the ligand

    """

    ec_df = pd.read_table(infile, delim_whitespace=True,
                          names=['pdb_template', 'tm_score', 'rmsd', 'seq_ident', 'seq_coverage',
                                 'c_score', 'ec_number', 'binding_residues'])
    ec_df['pdb_template_id'] = ec_df['pdb_template'].apply(lambda x: x[:4])
    ec_df['pdb_template_chain'] = ec_df['pdb_template'].apply(lambda x: x[4])

    ec_df = ec_df[['pdb_template_id', 'pdb_template_chain', 'tm_score', 'rmsd',
                   'seq_ident', 'seq_coverage', 'c_score', 'ec_number', 'binding_residues']]
    ec_df['c_score'] = pd.to_numeric(ec_df.c_score, errors='coerce')
    return ec_df  # .set_index('pdb_template_id')


def parse_go(infile):
    """Parse a GO output file from COACH and returns a dataframe

    The columns in all files are: GO terms, Confidence score, Name of GO terms.
        GO_MF.dat - GO terms in 'molecular function'
        GO_BP.dat - GO terms in 'biological process'
        GO_CC.dat - GO terms in 'cellular component'

    Args:
        infile (str): Path to any GO file

    Returns:
        go_df (Pandas DataFrame): Organized dataframe of results, columns defined below
            go_id: GO term ID
            go_term: GO term text
            c_score: confidence score of the GO prediction

    """

    go_df = pd.DataFrame()

    with open(infile) as go_file:
        for line in go_file.readlines():
            go_dict = {}

            go_split = line.split()
            go_dict['go_id'] = go_split[0]
            go_dict['c_score'] = go_split[1]
            go_dict['go_term'] = ' '.join(go_split[2:])

            go_df = go_df.append(go_dict, ignore_index=True)

    if len(go_df) == 0:
        return go_df

    go_df = go_df[['go_id', 'go_term', 'c_score']]
    go_df['c_score'] = pd.to_numeric(go_df.c_score, errors='coerce')
    return go_df  # .set_index('go_id')

class ITASSERParse():
    '''
    Parse I-TASSER & COACH results
    '''

    def __init__(self):
        pass
