#!/usr/bin/env python3

import collections
import pickle
from scipy.stats import hypergeom

### COMPARISON ###

def choose_valid_minimizers(minimizer_frequency, barcodes_minimizers, lim):
    #add data
    all_valid_mins = set()
    #end add data
    valid_minimizers = collections.defaultdict(list)
    for bc, mins in barcodes_minimizers.items():
        for m in set(mins):
            if minimizer_frequency[m] < lim:
                valid_minimizers[bc].append(m)
                #add data
                all_valid_mins.add(m)
                #end add data
    return valid_minimizers, all_valid_mins

def choose_valid_barcodes(valid_minimizers, lim):
    valid_barcodes = []
    for bc, mins in valid_minimizers.items():
        if len(mins) > lim:
            valid_barcodes.append(bc)
    return valid_barcodes

def set_footprints(valid_barcodes, valid_minimizers, lim_big, lim_small):
    bcm_big = {}
    bcm_small = {}
    for bc in valid_barcodes:
        mins = valid_minimizers[bc]
        if len(mins) > lim_big:
            m = sorted(mins)[:lim_big]
            bcm_big[bc] = set(m)
        else:
            m = sorted(mins)[:lim_small]
            bcm_small[bc] = set(m)
    return bcm_big, bcm_small

def compare_footprints(setA, setB, valid_minimizers_number):
    common_part = {}
    for bca, a in setA.items():
        for bcb, b in setB.items():
            cp = a.intersection(b)
            if cp:
                common_part[bca] = {}
                common_part[bca][bcb] = {}
                common_part[bca][bcb]["cp"] = cp
                common_part[bca][bcb]["pvalue"] = hypergeom.sf(len(cp), valid_minimizers_number, len(a), len(b))
    return common_part

with open('/Users/khaili/all_mins.pickle', 'rb') as f:
    all_mins = pickle.load(f)

with open('/Users/khaili/bc_mins.pickle', 'rb') as f:
    bc_mins = pickle.load(f)

#use bc_mins and all_mins from minimize.py
mins_freq = {m:all_mins.count(m) for m in all_mins}
valid_mins, vm_set = choose_valid_minimizers(mins_freq, bc_mins, 10000)

#how many minimizers are less frequent than 10000
vm_number = len(vm_set)
valid_bc = choose_valid_barcodes(valid_mins, 30)

bcm_big, bcm_small = set_footprints(valid_bc, valid_mins, 100, 30)
cp_dict = compare_footprints(bcm_small, bcm_big, vm_number)

print(cp_dict)
