import pysam
import dbm
from itertools import groupby
from biothings.utils.dataload import dict_sweep, unlist, value_convert, merge_duplicate_rows
from utils.hgvs import get_hgvs_from_vcf
# tabix file links from CADD http://cadd.gs.washington.edu/download

# number of fields/annotations
VALID_COLUMN_NO = 116
cadd_file_path = '/opt/myvariant.info/load_archive/cadd/whole_genome_SNVs_inclAnno.tsv.gz'

# convert one snp to json
def _map_line_to_json(fields):
    assert len(fields) == VALID_COLUMN_NO
    chrom = fields[0]
    chromStart = fields[1]
    ref = fields[2]
    alt = fields[4]
    HGVS = get_hgvs_from_vcf(chrom, chromStart, ref, alt)

    # load as json data
    if HGVS is None:
        return
    one_snp_json = {
        "_id": HGVS,
        "cadd": {
            'chrom': fields[0],
            'pos': fields[1],
            'ref': fields[2],
            'anc': fields[3],
            'alt': fields[4],
            'type': fields[5],
            'length': fields[6],
            'istv': fields[7],
            'isderived': fields[8],
            'annotype': fields[9],
            'consequence': fields[10],
            'consscore': fields[11],
            'consdetail': fields[12],
            'gc': fields[13],
            'cpg': fields[14],
            'mapability': {
                '20bp': fields[15],
                '35bp': fields[16]
            },
            'scoresegdup': fields[17],
            'phast_cons': {
                'primate': fields[18],
                'mammalian': fields[19],
                'vertebrate': fields[20]
            },
            'phylop': {
                'primate': fields[21],
                'mammalian': fields[22],
                'vertebrate': fields[23]
            },
            'gerp': {
                'n': fields[24],
                's': fields[25],
                'rs': fields[26],
                'rs_pval': fields[27]
            },
            'bstatistic': fields[28],
            'mutindex': fields[29],
            'dna': {
                'helt': fields[30],
                'mgw': fields[31],
                'prot': fields[32],
                'roll': fields[33]
            },
            'mirsvr': {
                'score': fields[34],
                'e': fields[35],
                'aln': fields[36]
            },
            'targetscans': fields[37],
            'fitcons': fields[38],
            'chmm': {
                'tssa': fields[39],
                'tssaflnk': fields[40],
                'txflnk': fields[41],
                'tx': fields[42],
                'txwk': fields[43],
                'enh': fields[44],
                # 'enh': fields[45],
                'znfrpts': fields[46],
                'het': fields[47],
                'tssbiv': fields[48],
                'bivflnk': fields[49],
                'enhbiv': fields[50],
                'reprpc': fields[51],
                'reprpcwk': fields[52],
                'quies': fields[53],
            },
            'encode': {
                'exp': fields[54],
                'h3k27ac': fields[55],
                'h3k4me1': fields[56],
                'h3k4me3': fields[57],
                'nucleo': fields[58],
                'occ': fields[59],
                'p_val': {
                    'comb': fields[60],
                    'dnas': fields[61],
                    'faire': fields[62],
                    'polii': fields[63],
                    'ctcf': fields[64],
                    'mycp': fields[65]
                },
                'sig': {
                    'dnase': fields[66],
                    'faire': fields[67],
                    'polii': fields[68],
                    'ctcf': fields[69],
                    'myc': fields[70]
                },
            },
            'segway': fields[71],
            'motif': {
                'toverlap': fields[72],
                'dist': fields[73],
                'ecount': fields[74],
                'ename': fields[75],
                'ehipos': fields[76],
                'escorechng': fields[77]
            },
            'tf': {
                'bs': fields[78],
                'bs_peaks': fields[79],
                'bs_peaks_max': fields[80]
            },
            'isknownvariant': fields[81],
            'esp': {
                'af': fields[82],
                'afr': fields[83],
                'eur': fields[84]
            },
            '1000g': {
                'af': fields[85],
                'asn': fields[86],
                'amr': fields[87],
                'afr': fields[88],
                'eur': fields[89]
            },
            'min_dist_tss': fields[90],
            'min_dist_tse': fields[91],
            'gene': {
                'gene_id': fields[92],
                'feature_id': fields[93],
                'ccds_id': fields[94],
                'genename': fields[95],
                'cds': {
                    'cdna_pos': fields[96],
                    'rel_cdna_pos': fields[97],
                    'cds_pos': fields[98],
                    'rel_cds_pos': fields[99]
                },
                'prot': {
                    'protpos': fields[100],
                    'rel_prot_pos': fields[101],
                    'domain': fields[102]
                }
            },
            'dst2splice': fields[103],
            'dst2spltype': fields[104],
            'exon': fields[105],
            'intron': fields[106],
            'oaa': fields[107],   # ref aa
            'naa': fields[108],   # alt aa
            'grantham': fields[109],
            'polyphen': {
                'cat': fields[110],
                'val': fields[111]
            },
            'sift': {
                'cat': fields[112],
                'val': fields[113]
            },
            'rawscore': fields[114],    # raw CADD score
            'phred': fields[115]        # log-percentile of raw CADD score
        }
    }

    obj = dict_sweep(unlist(value_convert(one_snp_json)), ["NA"])
    yield obj


def load_contig(contig):
    """contig is #chrm, from 1 to 23, X and Y"""
    tabix = pysam.Tabixfile(cadd_file_path)
    data = fetch_generator(tabix, contig)
    for doc in data:
        yield doc

def fetch_generator(tabix, contig):
    dbfile_path = 'home/kevinxin/cadd/' + 'cadd_id' + contig
    db = dbm.open(dbfile_path)
    ids = db.keys()
    set_ids = set(ids)
    print(len(ids))
    fetch = tabix.fetch(contig)
    rows = map(lambda x: x.split('\t'), fetch)
#   looking for annotype as 'codingtranscript', 'noncodingtranscript'
    annos = (row for row in rows if "CodingTranscript" in row[9] or
             get_hgvs_from_vcf(row[0], row[1], row[2], row[4]) in set_ids)
    json_rows = map(_map_line_to_json, annos)
    json_rows = (row for row in json_rows if row)
    row_groups = (it for (key, it) in groupby(json_rows, lambda row: row["_id"]))
    return (merge_duplicate_rows(rg, "cadd") for rg in row_groups)

def load_data(data_folder):
    for contig in [i for i in range(1,24)] + ["X","Y"]:
        load_contig(config)

