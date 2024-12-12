import csv
from biothings.utils.dataload import unlist
from biothings.utils.dataload import value_convert_to_number
from biothings.utils.dataload import merge_duplicate_rows
from utils.hgvs import get_hgvs_from_vcf
from itertools import groupby
from math import nan

VALID_COLUMN_NO = 31

'''this parser is for denovo-db v1.5 downloaded from
http://denovo-db.gs.washington.edu/denovo-db/Download.jsp'''


# convert one snp to json
def _map_line_to_json(df):
    # specific variable treatment
    chrom = df["Chr"]
    if chrom == 'M':
        chrom = 'MT'

    position = int(df["Position"])
    ref, alt = df["Variant"].upper().split(">")

    HGVS = get_hgvs_from_vcf(chrom, position, ref, alt, mutant_type=False)

    sampleid = df["SampleID"]
    studyname = df["StudyName"]
    pubmedid = df["PubmedID"]
    numprobands = df["NumProbands"]
    numcontrols = df["NumControls"]
    sequencetype = df["SequenceType"]
    primaryphenotype = df["PrimaryPhenotype"]
    validation = df["Validation"]
    chrom = df["Chr"]
    position = df["Position"]
    variant = df["Variant"]
    rsid = clean_rsid(df["rsID"], ("0", ))
    dbsnpbuild = clean_data(df["DbsnpBuild"], ("0", ))
    ancestralallele = df["AncestralAllele"]
    kgenomecount = df["1000GenomeCount"]
    exacfreq = df["ExacFreq"]
    espaafreq = df["EspAaFreq"]
    espeafreq = df["EspEaFreq"]
    transcript = clean_data(df["Transcript"], ("none", ""))
    codingdnasize = clean_data(df["codingDnaSize"], ("-1", ))
    gene = clean_data(df["Gene"], ("NA", ""))
    functionclass = clean_data(df["FunctionClass"], ("none", ""))
    cdnavariant = clean_data(df["cDnaVariant"], ("NA", ""))
    proteinvariant = clean_data(df["ProteinVariant"], ("NA", ""))
    exon_intron = clean_data(df["Exon_Intron"], ("NA",))
    polyphen_hdiv = clean_data(df["PolyPhen_HDiv"], ("-1",))
    polyphen_hvar = clean_data(df["PolyPhen_HVar"], ("-1",))
    siftscore = clean_data(df["SiftScore"], ("-1",))
    caddscore = clean_data(df["CaddScore"], ("-1",))
    lofscore = clean_data(df["LofScore"], ("-1",))
    lrtscore = clean_data(df["LrtScore"], ("-1",))

# load as json data
    one_snp_json = {
        "_id": HGVS,
        "denovodb": {
            "ref": ref,
            "alt": alt,
            "sampleid": sampleid,
            "studyname": studyname,
            "pubmedid": pubmedid,
            "numprobands": numprobands,
            "numcontrols":  numcontrols,
            "sequencetype":  sequencetype,
            "primaryphenotype": primaryphenotype,
            "validation": validation,
            "position": position,
            "variant": variant,
            "rsid": rsid,
            "dbsnpbuild": dbsnpbuild,
            "ancestralallele": ancestralallele,
            "1000genomecount":  kgenomecount,
            "exacfreq":  exacfreq,
            "espaafreq":   espaafreq,
            "espeafreq":   espeafreq,
            "transcript":  transcript,
            "codingdnasize":  codingdnasize,
            "gene":    gene,
            "functionclass":  functionclass,
            "cdnavariant": cdnavariant,
            "proteinvariant":  proteinvariant,
            "exon_intron":   exon_intron,
            "polyphen_hdiv": polyphen_hdiv,
            "polyphen_hvar": polyphen_hvar,
            "siftscore":   siftscore,
            "caddscore": caddscore,
            "lofscore":  lofscore,
            "lrtscore":  lrtscore,
        }
    }
    # one_snp_json = dict_sweep(unlist(value_convert_to_number(one_snp_json)), vals=[np.nan])

    one_snp_json = value_convert_to_number(one_snp_json)
    return one_snp_json


def clean_index(s):
    return s.replace("/", "_").replace("(", "_").replace(")", "").replace("#", "")


def clean_data(d, vals):
    if d in vals:
        return nan
    else:
        return d


def clean_rsid(d, vals):
    if d in vals:
        return nan
    else:
        return "rs{}".format(d)


# open file, parse, pass to json mapper
def load_data(input_file, version='hg19'):
    open_file = open(input_file)
    db_denovodb = csv.reader(open_file, delimiter="\t")
    index = next(db_denovodb)
    while index[0].startswith("##"):
        index = next(db_denovodb)
    assert len(index) == VALID_COLUMN_NO, "Expecting %s columns, but got %s" % (VALID_COLUMN_NO, len(index))
    index = [clean_index(s) for s in index]
    denovodb = (dict(zip(index, row)) for row in db_denovodb)
    denovodb = filter(lambda row: row["Chr"] != "", denovodb)
    json_rows = map(_map_line_to_json, denovodb)
    json_rows = (row for row in json_rows if row)
    json_rows = sorted(json_rows, key=lambda row: row["_id"])
    row_groups = (it for (key, it) in groupby(json_rows, lambda row: row["_id"]))
    json_rows = (merge_duplicate_rows(rg, "denovodb") for rg in row_groups)
    return (unlist(dict_sweep(row, vals=[nan, ])) for row in json_rows)
    # return (merge_duplicate_rows(rg, "denovodb") for rg in row_groups)


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
