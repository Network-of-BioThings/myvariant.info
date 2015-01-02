# -*- coding: utf-8 -*-
import csv
from itertools import groupby, imap, ifilter
from utils.dataload import dict_sweep, list_split, unlist, value_convert


VALID_COLUMN_NO = 70
    

# convert one snp to json
def _map_line_to_json(fields):
    assert len(fields) == VALID_COLUMN_NO
    chrom = fields[5]
    chromStart = fields[6]
    alleles = fields[58].split("/")
    allele1 = alleles[0]
    allele2 = alleles[1][0]
    HGVS = "chr%s:g.%s%s>%s" % (chrom, chromStart, allele1, allele2)

    # load as json data
    if HGVS is None:
        return

    one_snp_json = {

        "_id": HGVS,
        "grasp":
            {
                'hg19':
                     {
                         'chr': fields[5],
                         'pos': fields[6]
                     },
                 'nhlbi_key': fields[0],
                 'hupfield': fields[1],
                 'last_curation_date': fields[2],
                 'creation_date': fields[3],
                 'snpid_dbsnp134': fields[4],
                 'publication':
                     {
                         'pmid': fields[7],
                         'snpid_in_paper': fields[8],
                         'location_within_paper': fields[9],
                         'p_value': fields[10],
                         'phenotype': str(fields[11]),
                         'paper_phenotype_description': fields[12],
                         'paper_phenotype_categories': fields[13],
                         'date_pub': fields[14]
                     },
                 'in_nhgri_cat': fields[15],
                 'journal': fields[16],
                 'title': fields[17],
                 'includes_male_female_only_analyses': fields[18],
                 'exclusively_male_female': fields[19],
                 'initial_sample_description': fields[20],
                 'replication_sample_description': fields[21],
                 'platform_snps_passing_qc': fields[22],
                 'gwas_ancestry_description': fields[23],
                 'discovery':
                     {
                         'total_samples': fields[25],
                         'european': fields[26],
                         'african': fields[27],
                         'east_asian': fields[28],
                         'indian_south_asian': fields[29],
                         'hispanic': fields[30],
                         'native': fields[31],
                         'micronesian': fields[32],
                         'arab_me': fields[33],
                         'mixed': fields[34],
                         'unspecified': fields[35],
                         'filipino': fields[36],
                         'indonesian': fields[37]
                     },
                 'replication':
                     {
                         'total_samples': fields[38],
                         'european': fields[39],
                         'african': fields[40],
                         'east_asian': fields[41],
                         'indian_south_asian': fields[42],
                         'hispanic': fields[43],
                         'native': fields[44],
                         'micronesian': fields[45],
                         'arab_me': fields[46],
                         'mixed': fields[47],
                         'unspecified': fields[48],
                         'filipino': fields[49],
                         'indonesian': fields[50]
                     },
                 'in_gene': fields[51],
                 'nearest_gene': fields[52],
                 'in_lincrna': fields[53],
                 'in_mirna': fields[54],
                 'in_mirna_bs': fields[55],
                 'dbsnp':
                     {
                         'fxn': fields[56],
                         'maf': fields[57],
                         'alleles_het_se': fields[58],
                         'validation': fields[59],
                         'clin_status': fields[60],
                     },
                 'oreg_anno': fields[61],
                 'conserv_pred_tfbs': fields[62],
                 'human_enhancer': fields[63],
                 'rna_edit': fields[64],
                 'polyphen2': fields[65],
                 'sift': fields[66],
                 'ls_snp': fields[67],
                 'uniprot': fields[68],
                 'eqtl_meth_metab_study': fields[69]
            }
        }
    return list_split(dict_sweep(unlist(value_convert(one_snp_json)), [""]), ",")

## replace None indices with ''
def row_generator(db_row):
    ind = range(VALID_COLUMN_NO)
    row = []
    for i in ind:
        try:
            row.append(db_row[i])
        except:
            row.append('')
    return row
    
    
# open file, parse, pass to json mapper
def load_data(input_file):
    open_file = open(input_file)
    open_file = csv.reader(open_file, delimiter="\t")
    open_file.next()
    grasp = imap(row_generator, open_file)
    grasp = ifilter(lambda row: row[58] != "", grasp)
    json_rows = imap(_map_line_to_json, grasp)
    json_rows = (row for row in json_rows if row)
    row_groups = (it for (key, it) in groupby(json_rows, lambda row: row["_id"]))
    return (merge_duplicate_rows(rg, "grasp") for rg in row_groups)
