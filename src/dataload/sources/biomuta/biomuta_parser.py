import csv
from biothings.utils.dataload import unlist
from biothings.utils.dataload import value_convert_to_number
from biothings.utils.dataload import merge_duplicate_rows
from utils.hgvs import get_hgvs_from_vcf
from itertools import groupby

VALID_COLUMN_NO = 17

'''this parser is for BioMuta v3.0(BioMuta3 Complete Dataset) downloaded from
https://hive.biochemistry.gwu.edu/cgi-bin/prd/biomuta/servlet.cgi'''


# convert one snp to json
def _map_line_to_json(df):
    # specific variable treatment
    genomic_position = clean_data(df["genomic_position"],  ("-",))
    if not genomic_position:
        return
    genomic_position_split = genomic_position.replace("chr", "").replace("-", ":").split(":")
    chrom = genomic_position_split[0]
    chromStart = genomic_position_split[1]
    if chrom == 'M':
        chrom = 'MT'

    ref = df["ref_nuc"]
    alt = df["var_nuc"]

    HGVS = get_hgvs_from_vcf(chrom, int(chromStart), ref, alt, mutant_type=False)

    index = df["index"]
    uniprotkb_swiss_prot_id = clean_data(df["uniprotkb_swiss_prot_id"], ("-",))
    gene_name = clean_data(df["gene_name"], ("-",))
    refseq_nucleotide_id = clean_data(df["refseq_nucleotide_id"], ("-",))
    position_nuc = clean_data(df["position_nuc"], ("-",))
    position_aa = clean_data(df["position_aa"], ("-",))
    ref_aa = clean_data(df["ref_aa"], ("-",))
    var_aa = clean_data(df["var_aa"], ("-",))
    polyphen = clean_data(df["polyphen"], ("-",))
    pmid = clean_data(df["pmid"], ("-",))

    cancer_type = clean_data(df["cancer_type"], ("-",))
    if cancer_type:
        cancer_type_split = cancer_type.replace(" / ", ":").split(":")
        assert len(cancer_type_split) == 3, "cancer_type split error : {} : {}".format(HGVS, cancer_type)
        _d, doid, term = cancer_type_split
        assert _d == "DOID", "cancer_type split error : {} : {}".format(HGVS, cancer_type)
    else:
        doid = None
        term = None

    source = clean_data(df["source"], ("-",))
    vfunction = clean_data(df["function"], ("-",))
    if vfunction:
        vfunction = vfunction.split("|")
    status = clean_data(df["status"], ("-",))

# load as json data
    one_snp_json = {
        "_id": HGVS,
        "biomuta": {
            'index': index,
            'uniprotkb_swiss_prot_id': uniprotkb_swiss_prot_id,
            'gene_name': gene_name,
            'refseq_nucleotide_id': refseq_nucleotide_id,
            'genomic_position': genomic_position,
            'position_nuc': position_nuc,
            'ref_nuc': ref,
            'var_nuc': alt,
            'position_aa': position_aa,
            'ref_aa': ref_aa,
            'var_aa': var_aa,
            'polyphen': polyphen,
            'pmid': pmid,
            'cancer_type': {
                "DOID": doid,
                "term": term },
            'source': source,
            'function': vfunction,
            'status': status,
        }
    }
    one_snp_json = value_convert_to_number(one_snp_json)
    return one_snp_json


def clean_index(s):
    return s.lower().replace("/", "_").replace("-", "_").replace("(", "_").replace(")", "").replace("#", "")


def clean_data(d, vals):
    if d in vals:
        return None
    else:
        return d

# open file, parse, pass to json mapper
def load_data(input_file, version='hg19'):
    open_file = open(input_file)
    db_biomuta = csv.reader(open_file)
    index = next(db_biomuta)
    assert len(index) == VALID_COLUMN_NO, "Expecting %s columns, but got %s" % (VALID_COLUMN_NO, len(index))
    index = [clean_index(s) for s in index]
    biomuta = (dict(zip(index, row)) for row in db_biomuta)
    biomuta = filter(lambda row: row["index"] != "", biomuta)
    json_rows = map(_map_line_to_json, biomuta)
    json_rows = (row for row in json_rows if row)
    json_rows = sorted(json_rows, key=lambda row: row["_id"])
    row_groups = (it for (key, it) in groupby(json_rows, lambda row: row["_id"]))
    json_rows = (merge_duplicate_rows(rg, "biomuta") for rg in row_groups)
    return (unlist(dict_sweep(row, vals=[None, ])) for row in json_rows)


def dict_sweep(d, vals=[".", "-", "", "NA", "none", " ", "Not Available", "unknown"]):
    """
    @param d: a dictionary
    @param vals: a string or list of strings to sweep
    """
    for key, val in list(d.items()):
        if val in vals:
            del d[key]
        elif isinstance(val, list):
            val = [v for v in val if v not in vals]
            for item in val:
                if isinstance(item, dict):
                    dict_sweep(item, vals)
            if len(val) == 0:
                del d[key]
            else:
                d[key] = val
        elif isinstance(val, dict):
            dict_sweep(val, vals)
            if len(val) == 0:
                del d[key]
    return d