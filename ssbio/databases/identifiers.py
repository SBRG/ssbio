from bioservices.uniprot import UniProt
import requests

bsup = UniProt()


def kegg_mapper(kegg_organism, map_db='uniprot'):
    '''
    Generates a dictionary that maps KEGG gene IDs to NCBI Entrez gene IDs ('ncbi-geneid') or UniProt ('uniprot').
    Input:  kegg_organism - the KEGG organism name which you can determine from http://www.genome.jp/kegg/catalog/org_list.html
            map_db - ncbi-geneid OR uniprot (default): the database you want to map to
    Output: id_mapper - a dictionary of {KEGG_ID: mapped_ID}
    '''

    id_mapper = {}
    r = requests.post('http://rest.kegg.jp/conv/%s/%s' % (map_db,kegg_organism))

    for line in r.content.split('\n'):
        if not line: continue

        idents = line.split('\t')

        orig = idents[0].split(':')
        # orig_database = orig[0]
        orig_id = orig[1]

        conv = idents[1].split(':')
        # conv_database = conv[0]
        conv_id = conv[1]

        id_mapper[orig_id] = conv_id

    return id_mapper


def bioservices_uniprot_mapper(map_from, map_to, ident):
    return dict(bsup.mapping(fr=map_from, to=map_to, query=ident))