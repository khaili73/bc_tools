#!/usr/bin/env python3

import gzip
import collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import numpy as np
import argparse
from contextlib import ExitStack
import re
import os
import pickle

#fastq.gz with barcodes location: ~/data/barcodes1/outs/barcoded.fastq.gz 1-8
parser = argparse.ArgumentParser(description='Processing barcodes from fastq records', epilog="Used libraries: gzip, collections, Bio, numpy, argparse.")
###for local testing###
#parser.add_argument('infiles', metavar=('input'), help='input fastq.gz files with reads', nargs='+')
parser.add_argument('barcodes_dir', metavar=('bc_dir'), help='path to directory with longranger results')
parser.add_argument('barcoded_files', metavar=('bc_files'), help='relative path with barcoded files in "barcodesx" directory')
###END###
parser.add_argument('--x', default=8, help='number of directories with barcodes output')
parser.add_argument('--out', default=os.environ['HOME'], help='output files directory')
#parser.add_argument('--mins', default='/minimizers.tsv', help='output minimizers file name')
#parser.add_argument('--mbc', default='/mins_per_bc.tsv', help='output minimizers per barcode count file name')
parser.add_argument('--modulo', type=int, help='modulo value for minimizers creation')
args = parser.parse_args()

### MINIMIZERS ###

def get_minimizer(k, sequence):
    sequence.alphabet = generic_dna
    rc_sequence = sequence.reverse_complement
    minimizer = abs(hash(str(sequence)[:3]))
    for seq in str(sequence), str(rc_sequence):
        for i in range(len(seq) - k + 1):
            hashed_kmer = hash(seq[i:i+k])
            if args.modulo:
                hashed_kmer = hashed_kmer % args.modulo #negative values
            if abs(hashed_kmer) <= minimizer: #choose rigthmost lowest value
                minimizer = abs(hashed_kmer)
                rem_i = i
                rem_seq = seq
    return minimizer

def get_barcode(record):
    if re.search("BX:Z:[ATCG]{16}-1", record.description):
        return re.search("BX:Z:[ATCG]{16}-1", record.description).group()
    else:
        return "BX:Z:0"

def process_barcode(current_barcode_records):
    minimizers = []
    for rec in current_barcode_records:
        minimizers.append(get_minimizer(21, rec.seq))
    minimizers.sort()
    return minimizers
    #return np.array(minimizers, dtype="uint16")

def next_valid(parser):
    b = "BX:Z:0"
    while b == "BX:Z:0":
        current_record = next(parser)
        b = get_barcode(current_record)
    return current_record

def clean_records(read_parsers):
    new_parsers = []
    current_records = []
    for parser in read_parsers:
        try:
            current_records.append(next_valid(parser))
            new_parsers.append(parser)
        except StopIteration:
            pass
    return new_parsers, current_records

def parse_barcodes(current_barcode, read_parsers, current_records):
    current_barcode_records = []
    new_records = []
    new_parsers = []
    for parser, record in zip(read_parsers, current_records):
        current_record = record
        try:
            while get_barcode(current_record) == current_barcode:
                current_barcode_records.append(current_record)
                current_record = next_valid(parser)
            new_records.append(current_record)
            new_parsers.append(parser)
        except StopIteration:
            pass
    return new_records, new_parsers, current_barcode_records

def create_barcodes_minimizers_map(read_parsers, current_records):
    bc_mins = {}
    #additional data
    all_mins = []
    #end additional data
    while read_parsers:
        current_barcode = min(get_barcode(record) for record in current_records)
        current_records, read_parsers, current_barcode_records = parse_barcodes(current_barcode, read_parsers, current_records)
        mins = process_barcode(current_barcode_records)
        bc_mins[current_barcode] = list(set(mins)) #unique mins
        #additional data
        all_mins.extend(list(set(mins))) #unique mins
        #end add data
    return bc_mins, all_mins

with ExitStack() as stack:
    ###for local testing###
    paths = [args.barcodes_dir + "/barcodes" + str(x) + args.barcoded_files for x in range(1,args.x+1)]
    files = [stack.enter_context(gzip.open(fname, "rt")) for fname in paths]
    #files = [stack.enter_context(gzip.open(fname, "rt")) for fname in args.infiles]
    ###END###
    read_parsers = [SeqIO.parse(f, "fastq") for f in files]
    read_parsers, current_records = clean_records(read_parsers)
    bc_mins, all_mins = create_barcodes_minimizers_map(read_parsers, current_records)

with open(args.out+"/all_mins.pickle", "wb") as f:
    pickle.dump(all_mins, f)

with open(args.out+"/bc_mins.pickle", "wb") as f:
    pickle.dump(bc_mins, f)

