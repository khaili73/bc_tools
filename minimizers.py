#!/Users/hekate/anaconda3/bin/python3.7

#import hashlib

dna = "ACGTACTTGTACACCGTACAAGGGCTAGAGGGCTAGGCTACTAAATCAG"

def get_kmers(k, dna):
   mers = []
   for i in range(len(dna) - k + 1):
      mers.append(dna[i:i+k])
   return mers

def get_minimizer(mers):
   h_mers = []
   for mer in mers:
      #h_mers.append(hashlib.md5(mer.encode()))
      h_mers.append(hash(mer))
   return min(h_mers)

mers = get_kmers(3,dna)
print(get_minimizer(mers))

