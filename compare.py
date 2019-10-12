#!/usr/bin/env python3

import pickle
from scipy.stats import hypergeom
from collections import Counter, namedtuple

all_mins_file = '/home/kwegrzyn/data/mod24_pickles/all_mins.pickle'
bc_mins_file = '/home/kwegrzyn/data/mod24_pickles/bc_mins.pickle'

### COMPARISON ###

def set_footprints(barcodes, barcodes_minimizers, number):
    footprints = {}
    for bc in barcodes:
        barcodes_minimizers[bc].sort(reverse=True)
        ftp = barcodes_minimizers[bc][:number]
        footprints[bc] = set(ftp)
    return footprints


def compare_footprints(mins_a, mins_b):
    common_part = mins_a.intersection(mins_b)
    if common_part:
        return common_part


def count_pvalue(mins_a, mins_b, common_part, valid_minimizers_number):
    return hypergeom.sf(len(common_part), valid_minimizers_number, len(mins_a), len(mins_b))


def undo_pickle(path):
    with open(path, 'rb') as f:
        obj = pickle.load(f)
    return obj


def do_pickle(path, obj):
    with open(path, "wb") as f:
        pickle.dump(obj, f)


def choose_frequent(minimizers, cutoff):
    minimizers_frequency = Counter(minimizers)
    return [min for min, frequency in minimizers_frequency.items() if frequency > cutoff]


def sift_minimizers(bcs_with_mins, black_list):
    sagnificant_barcodes = {}
    for bc, mins in bcs_with_mins.items():
        sagnificant_mins = [min for min in mins if min not in black_list]
        sagnificant_barcodes[bc] = sagnificant_mins
    return sagnificant_barcodes


def group_barcodes(barcodes_minimizers, s_sup, m_sup):
    small = []
    medium = []
    large = []
    group = namedtuple('group', ['small', 'medium', 'large'])
    for bc, mins in barcodes_minimizers.items():
        if len(mins) < s_sup:
            small.append(bc)
        elif len(mins) < m_sup:
            medium.append(bc)
        else:
            large.append(bc)
    return group(small, medium, large)


def list_to_str(alist):
    return ','.join([str(a) for a in alist])

def save_comparison(f_name, footprints_a, footprints_b, valid_minimizers_number):
    with open(f_name, 'w') as f:
        for a_bc, a_mins in footprints_a.items():
            for b_bc, b_mins in footprints_b.items():
                common_part = compare_footprints(a_mins, b_mins)
                if common_part:
                    pvalue = count_pvalue(a_mins, b_mins, common_part, valid_minimizers_number)
                    s_pvalue = str(round(pvalue, 4))
                    s_cp = list_to_str(common_part)
                    f.write("\t".join([a_bc, b_bc, s_pvalue, s_cp, "\n"]))
    return


def save_footprint_info(f_name, footprint, sagnificant_bcmins):
    with open(f_name, 'w') as f:
        for bc, ftp in footprint.items():
            all_mins = sagnificant_bcmins[bc]
            s_ftp = list_to_str(ftp)
            s_mins = list_to_str(all_mins)
            f.write("\t".join([bc, s_ftp, s_mins, "\n"]))
    return


mins = undo_pickle(all_mins_file)
mins_black_list = choose_frequent(mins, 1000)

do_pickle('/home/kwegrzyn/data/black_list.pickle', mins_black_list)

valid_minimizers_number = len(mins) - len(mins_black_list)

barcode_minimizers = undo_pickle(bc_mins_file)
sagnificant_bcmins = sift_minimizers(barcode_minimizers, mins_black_list)

group = group_barcodes(sagnificant_bcmins, 30, 100)

footprints_small = set_footprints(group.small, sagnificant_bcmins, -1)
footprints_medium = set_footprints(group.medium, sagnificant_bcmins, 30)
footprints_large = set_footprints(group.large, sagnificant_bcmins, 100)

ftps = (footprints_small, footprints_medium, footprints_large)
fnames = ('/home/kwegrzyn/data/footprint_S.tsv', '/home/kwegrzyn/data/footprint_M.tsv', '/home/kwegrzyn/data/footprints_L.tsv')
for ftp, fname in zip(ftps, fnames):
    save_footprint_info(fname, ftp, sagnificant_bcmins)

save_comparison('/home/kwegrzyn/data/comparison_ML.tsv', footprints_medium, footprints_large, valid_minimizers_number)
save_comparison('/home/kwegrzyn/data/comparison_SM.tsv', footprints_medium, footprints_small, valid_minimizers_number)
