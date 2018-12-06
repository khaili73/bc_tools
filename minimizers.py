#!/Users/hekate/anaconda3/bin/python3.6

import numpy as np

r = "ACGTACTTGTACACCGTACAAGGGCTAGAGGGCTAGGCTACTAAATCAG"
bc = "ACGTACT"

def get_minimizer(k, read):
   minimizer = 31
   for i in range(len(read) - k + 1):
      hashed_kmer = hash(read[i:i+k]) % 31 #jak dobrac odpowiednie modulo; co zrobic z wartosciami ujemnymi
      if hashed_kmer <= minimizer:
         minimizer = hashed_kmer
   return minimizer

reads = [(r, bc)]
minimizers = {}

for read, barcode in reads:
   if barcode not in minimizers.keys():
      minimizers[barcode] = [] #zmienic na array
   minimizers[barcode].append(get_minimizer(5,r))

   
print(minimizers)
