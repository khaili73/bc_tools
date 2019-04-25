#!/usr/bin/env python2.7

import matplotlib
matplotlib.use("AGG")

import matplotlib.pyplot as plt
import numpy as np


with open("./minimizers_mod216.tsv", "r") as f:
   minimizers = []
   lines = f.readlines()
   for line in lines:
      minimizers.append(line.strip())

with open("./mins_per_bc_mod216_uniq.tsv", "r") as ff:
   min_per_bc = []
   lines = ff.readlines()
   for line in lines:
      min_per_bc.append(line.strip().split()[1])

#max_m = int(max(minimizers))
#max_mbc = int(max(min_per_bc))
max_m = 1000
max_mbc = 100

#np.array(['1','2','3']).astype(np.float)
m = np.array(minimizers).astype(np.int)
mbc = np.array(min_per_bc).astype(np.int)

plt.figure(1)
#plt.hist(m)
n, bins, p = plt.hist(m, bins = range(0, max_m+int(max_m/10), int(max_m/10)))
plt.xticks(bins)
plt.title("Minimizers")
plt.xlabel("Minimizer value")
plt.ylabel("Frequency")

plt.savefig("m.png")

plt.figure(2)
#plt.hist(mbc)
n, bins, p = plt.hist(mbc, bins = range(0, max_mbc+1, int(max_mbc/10)))
locs, labels = plt.xticks()
plt.xticks(bins)
plt.title("Minimizers per barcode")
plt.xlabel("Number of minimizers")
plt.ylabel("Frequency")

plt.savefig("mperb.png")
