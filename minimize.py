#!/Users/hekate/anaconda3/bin/python3.6

import gzip
import collections
from Bio import SeqIO

barcodes_file = "barcodes.fastq.gz"

def get_minimizer(k, read):
    minimizer = 2**16
    for i in range(len(read) - k + 1):
       hashed_kmer = hash(read[i:i+k]) % 2**16 #problem with negative values
       if hashed_kmer <= minimizer: #rigthmost lowest value choosen
          minimizer = hashed_kmer
    return minimizer

with gzip.open(barcodes_file, "rt") as handle:
   barcodes_dict = collections.defaultdict(list) #dictionary of sequences mapped to their barcodes {bc:[seq1,seq2]}
   for record in SeqIO.parse(handle, "fastq"):
      barcode = record.description.split(' ')[1].split(":")[2] #line with reads description
      minimizer = get_minimizer(20, str(record.seq))
      barcodes_dict[barcode].append(minimizer)
print(barcodes_dict)
