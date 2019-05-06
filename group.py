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

#fastq.gz with barcodes location: ~/data/barcodes1/outs/barcoded.fastq.gz 1-8
parser = argparse.ArgumentParser(description='Processing barcodes from fastq records', epilog="Used libraries: gzip, collections, Bio, numpy, argparse.")
###for local testing###
#parser.add_argument('infiles', metavar=('input'), help='input fastq.gz files with reads', nargs='+')
parser.add_argument('barcodes_dir', metavar=('bc_dir'), help='path to directory with longranger results')
parser.add_argument('barcoded_files', metavar=('bc_files'), help='relative path with barcoded files in "barcodesx" directory')
###END###
parser.add_argument('--x', default=8, help='number of directories with barcodes output')
parser.add_argument('--out', default=os.environ['HOME'], help='output files directory')
parser.add_argument('--mins', default='/minimizers.tsv', help='output minimizers file name')
parser.add_argument('--mbc', default='/mins_per_bc.tsv', help='output minimizers per barcode count file name')
parser.add_argument('--modulo', type=int, default=2**16, help='modulo value for minimizers creation (default: 2**16)')
args = parser.parse_args()

def get_minimizer(k, read):
    minimizer = args.modulo
    dna = Seq(read, generic_dna)
    cdna = dna.reverse_complement()
    for seq in str(dna), str(cdna):
        last = minimizer
        for i in range(len(seq) - k + 1):
            hashed_kmer = hash(seq[i:i+k]) % args.modulo #problem with negative values
            if hashed_kmer <= minimizer: #rigthmost lowest value chosen
                minimizer = hashed_kmer
        if last < minimizer:
            minimizer = last
    return minimizer

def get_barcode(record):
    if re.search("BX:Z:[ATCG]{16}-1", record.description):
        return re.search("BX:Z:[ATCG]{16}-1", record.description).group()
    else:
        return "BX:Z:0"

def process_barcode(current_barcode_records):
    minimizers = []
    for rec in current_barcode_records:
        minimizers.append(get_minimizer(21, str(rec.seq)))
    minimizers.sort()
    return minimizers
    #return np.array(minimizers, dtype="uint16")

def next_valid(parser):
    b = "BX:Z:0"
    while b == "BX:Z:0":
        current_record = next(parser)
        b = get_barcode(current_record)
    return current_record

with ExitStack() as stack:
    ###for local testing###
    paths = [args.barcodes_dir + "/barcodes" + str(x) + args.barcoded_files for x in range(1,args.x+1)]
    files = [stack.enter_context(gzip.open(fname, "rt")) for fname in paths]
    #files = [stack.enter_context(gzip.open(fname, "rt")) for fname in args.infiles]
    ###END###
    bc_min_dict = collections.defaultdict(list)
    read_parsers = [SeqIO.parse(f, "fastq") for f in files]
    current_records = []
    new_parsers = []
    for parser in read_parsers:
        try:
            current_records.append(next_valid(parser))
            new_parsers.append(parser)
        except StopIteration:
            pass
    read_parsers = new_parsers
    while read_parsers:
        current_barcode = min(get_barcode(record) for record in current_records)
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
        current_records = new_records
        read_parsers = new_parsers
        with open("%s.fastq" % current_barcode, "w") as output_handle:
            SeqIO.write(current_barcode_records, output_handle, "fastq")
