#!/usr/bin/env python3

import gzip
import collections
from Bio import SeqIO
import numpy as np
import argparse
from contextlib import ExitStack
import re
from datetime import datetime
startTime = datetime.now()

#fastq.gz with barcodes location: ~/data/barcodes1/outs/barcoded.fastq.gz 1-8
parser = argparse.ArgumentParser(description='Processing barcodes from fastq records', epilog="Used libraries: gzip, collections, Bio, numpy, argparse.")
#parser.add_argument('infiles', metavar=('input'), help='input fastq.gz files with reads', nargs='+')
parser.add_argument('barcodes_dir', metavar=('bc_dir'), help='path to directory with longranger results')
parser.add_argument('barcoded_files', metavar=('bc_files'), help='relative path with barcoded files in "barcodesx" directory')
parser.add_argument('--modulo', type=int, default=2**16, help='modulo value for minimizers creation (default: 2**16)')
args = parser.parse_args()

def get_minimizer(k, read):
    minimizer = 2**16
    for i in range(len(read) - k + 1):
        hashed_kmer = hash(read[i:i+k]) % args.modulo #problem with negative values
        if hashed_kmer <= minimizer: #rigthmost lowest value chosen
            minimizer = hashed_kmer
    return minimizer

def barcode(record):
    if re.search("BX:Z:[ATCG]{16}-1", record.description):
        return re.search("BX:Z:[ATCG]{16}-1", record.description).group()
    else:
        return "BX:Z:0"

def process_barcode(current_barcode_records):
    minimizers = []
    for rec in current_barcode_records:
        minimizers.append(get_minimizer(21, str(rec.seq)))
    minimizers.sort()
    return np.array(minimizers, dtype="uint16")

#not used yet
def next_valid(parser):
    try:
        next(parser)
    except StopIteration:
        read_parsers.remove(parser)
    return

with ExitStack() as stack:
    paths = [args.barcodes_dir + "barcodes" + str(x) + args.barcoded_files for x in range(1,9)]
    files = [stack.enter_context(gzip.open(fname, "rt")) for fname in paths]
    print(files)
    bc_min_dict = collections.defaultdict(list)
    read_parsers = [SeqIO.parse(f, "fastq") for f in files]
    current_records = [next(parser) for parser in read_parsers]
    with open("minimizers.tsv", "w") as out_mins, open("mins_per_bc.tsv", "w") as mbc:
        while read_parsers:
            current_barcode = min(barcode(record) for record in current_records)
            current_barcode_records = []
            new_records = []
            for parser, record in zip(read_parsers, current_records):
                current_record = record
                while barcode(current_record) == current_barcode:
                    current_barcode_records.append(current_record)
                    try:
                        current_record=next(parser)
                    except StopIteration:
                        read_parsers.remove(parser)
                        break
                new_records.append(current_record)
            current_records = new_records
            mins = process_barcode(current_barcode_records)
            if current_barcode != "BX:Z:0":
                for m in mins:
                    out_mins.write(str(m)+"\n")
                mbc.write(current_barcode + "\t" + str(len(mins)) + "\n")

with open("summary.txt", "w+") as s:
    s.write("Programm terminated succesfully.")
    s.write("Time lapsed since begining: %s" % str(datetime.now() - startTime))

