#!/usr/bin/env python3.6

import gzip
import collections
from Bio import SeqIO
import numpy as np
import argparse
from contextlib import ExitStack

#fastq.gz z barkodami sa tu: ~/data/barcodes1/outs/barcoded.fastq.gz 1-8
parser = argparse.ArgumentParser(description='Processing barcodes from fastq records', epilog="Used libraries: gzip, collections, Bio, numpy, argparse.")
parser.add_argument('infiles', metavar=('input'), help='input fastq.gz files with reads', nargs='+')
args = parser.parse_args()

'''
numpy uint16 - 2 bajty w array
2^8 - 2^16 histogramy dla różnych modulo
'''
def get_minimizer(k, read):
    minimizer = 2**16
    for i in range(len(read) - k + 1):
       hashed_kmer = hash(read[i:i+k]) % 2**16 #problem with negative values
       if hashed_kmer <= minimizer: #rigthmost lowest value choosen
          minimizer = hashed_kmer
    return minimizer

def barcode(record):
   if "-1" in record.description:
      return record.description.split(' ')[1].split(':')[2] 
   else:
      return False

#def process_barcode(current_barcode, current_barcode_records, bc_min_dict):
def process_barcode(current_barcode_records):
   minimizers = []
   for rec in current_barcode_records:
      minimizers.append(get_minimizer(21, str(rec.seq)))
      #bc_min_dict[current_barcode].append(get_minimizer(20, str(rec.seq)))
   minimizers.sort()
   #bc_min_dict[current_barcode] = np.array(minimizers, dtype="uint16")
   #min_per_bc[current_barcode] = len(minimizers)
   #mini_v.extend(minimizers)
   #return bc_min_dict
   return np.array(minimizers, dtype="uint16")

with ExitStack() as stack:
   files = [stack.enter_context(gzip.open(fname, "rt")) for fname in args.infiles]
   bc_min_dict = collections.defaultdict(list) 
   read_parsers = [SeqIO.parse(f, "fastq") for f in files]
   current_records = [next(parser) for parser in read_parsers]
   with open("minimizers.tsv", "w") as out_mins, open("mins_per_bc.tsv", "w") as mbc:
      while read_parsers:
         valid_current_records = []
         for parser, record in zip(read_parsers, current_records):
            while barcode(record) == False:
               try:
                  record=next(parser)
               except:
                  read_parsers.remove(parser)
                  current_records.remove(record)
            if barcode(record):
               valid_current_records.append(record)
         current_records = valid_current_records
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
            #while barcode(current_record) == False:
            #   try:
            #      current_record=next(parser)
            #   except StopIteration:
            #      read_parsers.remove(parser)
            new_records.append(current_record)
         current_records = new_records
         #d = process_barcode(current_barcode, current_barcode_records, bc_min_dict)
         mins = process_barcode(current_barcode_records)
         for m in mins:
            out_mins.write(str(m)+"\n")
         mbc.write(current_barcode + "\t" + str(len(mins)) + "\n")

#generate plots

'''
f, ax = plt.subplots(2, 1)

ax[0].bar(min_per_bc.keys(), min_per_bc.values(), width = 0.2, color='g')
ax[0].set_xticks(bc_names)
ax[0].set_xticklabels(bc_names, rotation=65)
ax[1].hist(mini_v, color='b')

plt.show()
'''
