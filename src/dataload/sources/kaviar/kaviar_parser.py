import vcf
from biothings.utils.dataload import unlist
from biothings.utils.dataload import value_convert_to_number, merge_duplicate_rows, dict_sweep
from utils.hgvs import get_hgvs_from_vcf
from itertools import groupby, chain
import os, csv

VALID_COLUMN_NO = 8

'''this parser is for Kaviar version 160204-Public(All variants, annotated with data sources) downloaded from
http://db.systemsbiology.net/kaviar/Kaviar.downloads.html'''


# convert one snp to json
def _map_line_to_json(item):
    chrom = item.CHROM
    chromStart = item.POS
    ref = item.REF
    info = item.INFO

    try:
        af = info['AF']
    except:
        af = None
    try:
        ac = info['AC']
    except:
        ac = None
    try:
        an = info['AN']
    except:
        ac = None
    try:
        ds = info['DS']
    except:
        ds = None

    # convert vcf object to string
    item.ALT = [str(alt) for alt in item.ALT]

    # if multiallelic, put all variants as a list in multi-allelic field
    hgvs_list = None
    if len(item.ALT) > 1:
        hgvs_list = []
        for alt in item.ALT:
            try:
                hgvs_list.append(get_hgvs_from_vcf(chrom, chromStart, ref, alt, mutant_type=False))
            except:
                hgvs_list.append(alt)

        assert len(item.ALT) == len(info['AC']), "Expecting length of item.ALT= length of info.AC, but not for %s" % (item)
        assert len(item.ALT) == len(info['AF']), "Expecting length of item.ALT= length of info.AF, but not for %s" % (item)
        if ds:
            if len(item.ALT) != len(info['DS']):
                ds_str = ",".join(info['DS'])
                ds_str = ds_str.replace("NA7022,18", "NA7022_18")
                ds_list = ds_str.split(",")
                info['DS'] = [d.replace("NA7022_18", "NA7022,18") for d in ds_list]
                assert len(item.ALT) ==len(info['DS']), "info.DS mismatch, %s: %s\n## DS: %s" % (item, info['DS'])

    for i, alt in enumerate(item.ALT):
        try:
            (HGVS, var_type) = get_hgvs_from_vcf(chrom, chromStart, ref, alt, mutant_type=True)
        except:
            continue

        if HGVS is None:
            return

        # load as json data
        one_snp_json = {
            "_id": HGVS,
            "kaviar": {
                "multi-allelic": hgvs_list,
                "ref": ref,
                "alt": alt,
                "af": info['AF'][i],
                "ac": info['AC'][i],
                "an": an,
                "ds": info['DS'][i].split("|") if ds else None,
            }
        }

        yield value_convert_to_number(one_snp_json)


# open file, parse, pass to json mapper
def load_data(input_file):
    vcf_reader = vcf.Reader(open(input_file, 'r'), strict_whitespace=True)
    json_rows = map(_map_line_to_json, vcf_reader)
    json_rows = chain.from_iterable(json_rows)

    if not os.path.exists("alldata.csv"):
        print("Writing data")
        with open("alldata.csv", "w") as f:
            dbwriter = csv.writer(f)
            for doc in json_rows:
                dbwriter.writerow([doc['_id'], str(doc)])
    if not os.path.exists("sorted.csv"):
        print("Start Sorting")
        import subprocess
        p = subprocess.Popen('sort alldata.csv > sorted.csv', shell=True)
        os.waitpid(p.pid, 0)
        print("Sorted")

    json_rows = csv.reader(open('sorted.csv'))
    json_rows = (eval(row[1]) for row in json_rows)
    row_groups = (it for (key, it) in groupby(json_rows, lambda row: row["_id"]))
    json_rows = (merge_duplicate_rows(rg, "kaviar") for rg in row_groups)
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