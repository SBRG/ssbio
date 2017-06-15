import requests
from bs4 import BeautifulSoup
import ssbio.utils

### Server

# http://opm.phar.umich.edu/server.php
#
# ### Definitions:
#
# Depth or hydrophobic thickness. This parameter indicates the calculated hydrophobic thickness
# (for TM proteins) or maximal penetration depth of protein atoms into the lipid hydrocarbon core
# (for peripheral/monotopic) proteins. The ± values for the depth and tilt angle show fluctuations of
# the corresponding parameters within 1 kcal/mol around the global minimum of transfer energy.
#
# Tilt angle is calculated between membrane normal (Z axis) and protein axis. The protein axis is
# calculated as the sum of TM secondary structure segment vectors (for TM proteins) or as the principal
# inertia axis (for peripheral proteins).
#
# Transfer energy of the protein from water to lipid bilayer. This energy roughly corresponds to the
# actual membrane binding energy for water-soluble peripheral proteins, unless: (a) some membrane-anchoring
# elements are missing or disordered in the crystal structure; (b) there is strong specific binding of lipids
# (e.g. in PH domains and other "lipid clamps"), or (c) membrane binding is coupled with significant structural
# changes of the protein (e.g. helix-coil transition for amphiphilic α-helices). In situations (a) and (b), the
# calculated membrane binding free energy is underestimated. In situation (c) it is usually overestimated.
#
# Table of membrane-embedded residues consists of two parts: (a) list of residues penetrating into the hydrocarbon
# core of the lipid bilayer for each subunit (tilt angles of individual subunits in this part are calculated based
# on the principal inertia axis), and (b) parts of transmembrane alpha-helices or beta-strands that are embedded
# into the hydrocarbon core (the entire secondary structures can be longer; tilt angles of individual subunits
# in this part are calculated as vector averages of TM secondary structure segment vectors).
#
# Output coordinate file for a protein positioned in the lipid bilayer. The origin of coordinates corresponds
# to the center of lipid bilayer. Z axis coincides with membrane normal; atoms with the positive sign of Z coordinate
# are arranged in the "outer" leaflet as defined by the user-specified topology. Positions of DUMMY atoms correspond
# to locations of lipid carbonyl groups.
#
# Diagnostic messages. The server produces diagnostic messages in a separate window. Please report these messages to
# developer if program aborts.

## SEE: https://structure.dynamic.ucsd.edu:9998/notebooks/projects_unsynced/sandbox/PPM_server_test.ipynb

def run_ppm_server(pdb_file, outfile, force_rerun=False):
    if ssbio.utils.force_rerun(outfile=outfile, flag=force_rerun):
        url = 'http://sunshine.phar.umich.edu/upload_file.php'
        files = {'userfile': open(pdb_file, 'rb')}
        r = requests.post(url, files=files)
        info = r.text

        # Save results in raw HTML format
        with open(outfile, 'w') as f:
            f.write(info)

    else:
        # Utilize existing saved results
        with open(outfile, 'r') as f:
            info = f.read()

    # Clean up the HTML stuff
    t = info.replace('\n', '')
    tt = t.replace('\r', '')
    ttt = tt.replace('\t', '')

    soup = BeautifulSoup(ttt, "lxml")
    # Find all tables in the HTML code
    tables = soup.find_all("table", attrs={"class": "data"})

    info_dict = {}

    # There are multiple tables with information
    table_index = 0
    for t in tables:
        data_index = 0

        # "row1" contains data
        for data in t.find_all('tr', attrs={"class": "row1"}):
            data_list = list(data.strings)

            if table_index == 0:
                info_dict['Depth/Hydrophobic Thickness'] = data_list[0]
                info_dict['deltaG_transfer'] = data_list[2]
                info_dict['Tilt Angle'] = data_list[3]

            if table_index == 1 and data_index == 0:
                info_dict['Embedded_residues_Tilt'] = data_list[0]
                info_dict['Embedded_residues'] = data_list[1]

            if table_index == 1 and data_index == 1:
                info_dict['Transmembrane_secondary_structure_segments_Tilt'] = data_list[0]
                info_dict['Transmembrane_secondary_structure_segments'] = data_list[1]

            if table_index == 2:
                info_dict['Output Messages'] = data_list[1]

            if table_index == 3:
                baseurl = 'http://sunshine.phar.umich.edu/'
                a = data.find('a', href=True)
                download_url = baseurl + a['href'].replace('./', '')
                info_dict['Output file download link'] = download_url

            data_index += 1
        table_index += 1

    return info_dict