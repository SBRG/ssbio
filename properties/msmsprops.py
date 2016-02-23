from Bio import PDB
from loader import Loader
l = Loader()

def msms_output(filename):
    my_structure = l.structure_reader(filename)
    model = my_structure[0]
    rd = PDB.ResidueDepth(model, filename)

    akeys = list(rd)
    if len(akeys) == 0:
        akeys = [('<Residue UNKNOWN het=  resseq=0 icode= >', (0, 0))]

    anslist = []
    if akeys == [('NA', ('NA', 'NA'))]:
        anslist = ["NA", "NA", "NA", "NA"]
        # print "warning: error at index:"+str(len(msmslist))

    else:
        for j in akeys:
            chain = [x.id for x in PDB.Selection.unfold_entities(j[0], 'C')][0]
            seqnum = j[0].id[1]
            re_depth = j[1][0]
            ca_depth = j[1][1]
            if chain == ' ':
                chain = 'X'
            anslist.append([chain, seqnum, re_depth, ca_depth])

    return anslist

def residue_depth(anslist):
    '''Computes the average residue and CA depth
    Returns a dictionary of "redepth" and "cadepth" (floats)
    '''
    depth_info = {}

    if anslist != ['NA', 'NA', 'NA', 'NA']:
        redepth = [x[2] for x in anslist]
        cadepth = [x[3] for x in anslist]
        depth_info['redepth'] = sum(redepth) / len(redepth)
        depth_info['cadepth'] = sum(cadepth) / len(cadepth)
    else:
        depth_info['redepth'] = None
        depth_info['cadepth'] = None

    return depth_info

if __name__ == '__main__':
    import glob
    files = glob.glob('test_structures/*')
    for f in files:
        print(f)
        msms = msms_output(f)
        print(residue_depth(msms))
