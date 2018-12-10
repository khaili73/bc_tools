#!/Users/hekate/anaconda3/bin/python3.6

import numpy as np
import gzip
#import sys

def get_minimizer(k, read):
   minimizer = 2**16
   for i in range(len(read) - k + 1):
      hashed_kmer = hash(read[i:i+k]) % 2**16 #co zrobic z wartosciamiprint(read[i:i+k])
      if hashed_kmer <= minimizer: #ostatnia najmniejsza wartość
         minimizer = hashed_kmer
   return minimizer

minimizers = {}

with gzip.open('barcodes.fastq.gz') as bc:
   lines = 0
   old_barcode = None
   for line in bc:
      lines += 1
      if "BX:Z:".encode("utf-8") in line:
         has_barcode = True
         barcode = line.strip().split(' '.encode("utf-8"))[-1].split(":".encode("utf-8"))[-1].split("-".encode("utf-8"))[0]
         if barcode not in minimizers.keys():
            if old_barcode: #problem z pierwszym i ostatnim
               minimizers[old_barcode] = np.asarray(minimizers[old_barcode]).sort()
            minimizers[barcode] = []
      elif lines % 4 == 1:
         has_barcode = False
      if has_barcode and lines % 4 == 2:
         minimizers[barcode].append(get_minimizer(15,line.strip()))
      old_barcode = barcode

print(minimizers)

