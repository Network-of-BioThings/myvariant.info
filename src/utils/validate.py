from __future__ import print_function
import re
import os.path
import time

from bitarray import bitarray

from biothings.utils.common import loadobj, is_str, open_anyfile, timesofar
from biothings.utils.mongo import get_src_db, doc_feeder

from warnings import warn


def nuc_to_bit(sequence):
    '''encode nucleotide into bit form'''
    code = {'Y': bitarray('0000'),
            'A': bitarray('0001'),
            'C': bitarray('0010'),
            'G': bitarray('0011'),
            'T': bitarray('0100'),
            'g': bitarray('0011'),
            't': bitarray('0100'),
            'a': bitarray('0001'),
            'c': bitarray('0010'),
            'N': bitarray('0101'),
            'M': bitarray('0110'),
            'R': bitarray('0111'),
            'W': bitarray('1000'),
            'K': bitarray('1001'),
            'n': bitarray('0101'),
            }
    #assert sequence in code, "%s" % repr(sequence)
    seq_bit = bitarray()
    try:
        seq_bit.encode(code, sequence)
    except ValueError:
        print("%s" % repr(sequence))
        raise
    return(seq_bit)


def bit_to_nuc(sequence):
    '''encode bit form to nucleotide'''
    if sequence == bitarray('0001'):
        nuc = 'A'
    elif sequence == bitarray('0010'):
        nuc = 'C'
    elif sequence == bitarray('0011'):
        nuc = 'G'
    elif sequence == bitarray('0100'):
        nuc = 'T'
    elif sequence == bitarray('0110'):
        nuc = 'M'
    elif sequence == bitarray('0101'):
        nuc = 'N'
    elif sequence == bitarray('0111'):
        nuc = 'R'
    elif sequence == bitarray('1000'):
        nuc = 'W'
    else:
        raise ValueError("Cannot decode input bits: %s" % repr(sequence))
    return nuc


def bit_to_nuc2(bits):
    '''a util function to convert a encoded bitarray back to
       nt sequence.
    '''
    code = {'A': bitarray('001'),
            'C': bitarray('010'),
            'G': bitarray('011'),
            'T': bitarray('100'),
            'N': bitarray('101'),
            'M': bitarray('110'),
            'R': bitarray('111')}
    return bits.decode(code)


def parse(str):
    '''parse variant name, print the variant name and
       return chromosome number, nucleotide position
       and nucleotide name
    '''
    pat = 'chr(\w+):g\.(\d+)(\w)\>(\w)'
    mat = re.match(pat, str)
    if mat:
        r = mat.groups()
        return (r[0], r[1], r[2])


def get_genome_in_bit(chr_fa_folder):
    ''' encode each chromosome fasta sequence into a bitarray,
        and store them in a dictionary with chr numbers as keys
        chr_fa_folder is the folder to put all gzipped fasta files:

        fasta files can be downloaded from NCBI FTP site:

        ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37.p13/Primary_Assembly/assembled_chromosomes/FASTA/
        chr<i>.fa.gz  (e.g. chr1.fa.gz)

    '''
    chr_bit_d = {}
    chr_range = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
    t0 = time.time()
    for i in chr_range:
        t1 = time.time()
        #file_name = 'hs_ref_GRCh37.p5_chr{}.fa.gz'.format(i)
        file_name = 'chr{}.fa.gz'.format(i)
        print("Loading {}...".format(file_name), end='')
        file_name = os.path.join(chr_fa_folder, file_name)
        with open_anyfile(file_name) as seq_f:
            seq_f.readline()   # skip header
            seq_bit = bitarray()
            for line in seq_f:
                line = line.rstrip('\n')
                line_bit = nuc_to_bit(line)
                seq_bit += line_bit
            chr_bit_d.update({i: seq_bit})
        print("done.[{}]".format(timesofar(t1)))
    print('='*20)
    print("Finished. [{}]".format(timesofar(t0)))

    return chr_bit_d


class VariantValidator:
    def __init__(self):
        self._chr_data = None

    def load_chr_data(self,genome_file):
        print("\tLoading chromosome data...", end='')
        self._chr_data = loadobj(genome_file)
        print("Done.")

    def validate_hgvs(self, hgvs_id, verbose=False):
        '''validate single hgvs variant name, return True/False,
           or None if input hgvs_id cannot be validated (could be
           wrong format, or ins/del type currently we don't validate.
        '''
        r = parse(hgvs_id)
        # set the range for r[0]
        chr_range = [str(i) for i in range(1, 23)] + ['X', 'Y', 'M', 'MT']
        chr_range = set(chr_range)

        if r and (str(r[0]) in chr_range):
            # get the chromosome sequence in bit form
            if self._chr_data is None:
                print()
                self.load_chr_data()
            if r[0] == 'M':
                chr = 'MT'
            else:
                chr = r[0]
            pos = int(r[1])
            nuc_hgvs = r[2]

            chr_bit = bitarray()
            chr_bit = self._chr_data[str(chr)]

            # get the nucleotide in chromsome sequence in bit form
            nuc_chr_bit = bitarray()
            nuc_chr_bit = chr_bit[pos*3-3:pos*3]
            nuc_chr = bit_to_nuc(nuc_chr_bit)

            # compare HGVS id with genome
            matched = nuc_hgvs == nuc_chr
            if verbose:
                if matched:
                    print('"{}":\t{}'.format(hgvs_id, matched))
                else:
                    print('"{}":\t{} (should be "{}")'.format(hgvs_id, matched, nuc_chr))
            return matched
        else:
            if verbose:
                print('"{}":\tNone(not tested).'.format(hgvs_id))
            return None

    def validate_many(self, hgvs_li, verbose=False, summary=True):
        '''validate multiple hgvs variant name'''
        out = []
        for hgvs_id in hgvs_li:
            out.append(self.validate_hgvs(hgvs_id, verbose=verbose))

        if summary:
            # print out counts
            print("# of VALID HGVS IDs:\t{0}".format(len([x for x in out if x is True])))
            print("# of INVALID HGVS IDs:\t{0}".format(len([x for x in out if x is False])))
            print("# of HGVS IDs skipped:\t {0}".format(len([x for x in out if x is None])))
        return out

    def validate_src(self, collection, return_false=False,
                     return_none=False, return_true=False, verbose=False, flag_invalid=False, generator=False):
        '''Validate hgvs ids from a src collection.'''

        return_dict = {
            False: return_false,
            True: return_true,
            None: return_none
        }

        # read in the collection from mongodb
        if is_str(collection):
            src = get_src_db()
            _coll = src[collection]
        else:
            _coll = collection
        cursor = doc_feeder(_coll, step=10000)

        out = {}
        print_only = not (return_false or return_none or return_true)
        if not print_only:
            # output dictionary, three keys: 'false','true','none'
            for k in return_dict:
                if return_dict[k]:
                    out[k] = []

        # initialize the count
        cnt_d = {True: 0, False: 0, None: 0}    # cnt_d
        # validate each item in the cursor
        for item in cursor:
            _id = item['_id']
            valid = self.validate_hgvs(_id, verbose=verbose)
            if valid == False and flag_invalid:
                collection.update({"_id": _id}, {'$set':{"unmatched_ref": "True"}})
            cnt_d[valid] += 1
            if return_dict[valid]:
                out[valid].append(_id)

        # print out counts
        print("\n# of VALID HGVS IDs:\t{0}".format(cnt_d[True]))
        print("# of INVALID HGVS IDs:\t{0}".format(cnt_d[False]))
        print("# of HGVS IDs skipped:\t {0}".format(cnt_d[None]))

        out['summary'] = cnt_d
        return out

    def validate_generator(self, generator, return_false=False,
                     return_none=False, return_true=False, verbose=False, flag_invalid=False):
        '''Validate hgvs ids from a src collection.'''

        warn(
            RuntimeWarning("Method validate_generator is not correctly implemented. "
                           "See https://github.com/biothings/myvariant.info/blob/master/CONTRIBUTE.md#variant_id-validation"),
            stacklevel=2,  # unwind 2 calls in the callstack
        )

        return_dict = {
            False: return_false,
            True: return_true,
            None: return_none
        }

        out = {}
        print_only = not (return_false or return_none or return_true)
        if not print_only:
            # output dictionary, three keys: 'false','true','none'
            for k in return_dict:
                if return_dict[k]:
                    out[k] = []

        # initialize the count
        cnt_d = {True: 0, False: 0, None: 0}    # cnt_d
        # validate each item in the cursor
        for item in generator:
            _id = item['_id']
            valid = self.validate_hgvs(_id, verbose=verbose)
            if valid == False and flag_invalid:
                collection.update({"_id": _id}, {'$set':{"unmatched_ref": "True"}})
            cnt_d[valid] += 1
            if return_dict[valid]:
                out[valid].append(_id)

        # print out counts
        print("\n# of VALID HGVS IDs:\t{0}".format(cnt_d[True]))
        print("# of INVALID HGVS IDs:\t{0}".format(cnt_d[False]))
        print("# of HGVS IDs skipped:\t {0}".format(cnt_d[None]))

        out['summary'] = cnt_d
        return out
