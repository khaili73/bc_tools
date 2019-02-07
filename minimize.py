#!/Users/hekate/anaconda3/bin/python3.6

import gzip
import collections
from Bio import SeqIO

barcodes_file = "barcoded.fastq.gz", "barcodes2.fastq.gz"
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

def process_barcode(current_barcode, current_barcode_records, bc_min_dict):
   for rec in current_barcode_records:
      bc_min_dict[current_barcode].append(get_minimizer(20, str(rec.seq)))
   return bc_min_dict

with gzip.open(barcodes_file[0], "rt") as file1, gzip.open(barcodes_file[1], "rt") as file2:
   bc_min_dict = collections.defaultdict(list) 
   read_parsers = [SeqIO.parse(f, "fastq") for f in (file1, file2)]
   current_records = [next(parser) for parser in read_parsers]
   while read_parsers:
      current_barcode = min(barcode(record) for record in current_records if barcode(record)) 
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
      d = process_barcode(current_barcode, current_barcode_records, bc_min_dict)
