#!/usr/bin/env python2.7

import matplotlib
matplotlib.use("AGG")

import matplotlib.pyplot as plt
import numpy as np


with open("/home/kwegrzyn/tools/minimizers.tsv", "r") as f:
   minimizers = []
   lines = f.readlines()
   for line in lines:
      minimizers.append(line.strip())    

with open("/home/kwegrzyn/tools/mins_per_bc.tsv", "r") as ff:
   min_per_bc = []
   lines = ff.readlines()
   for line in lines:
      min_per_bc.append(line.strip().split()[1])

max_m = int(max(minimizers))
max_mbc = int(max(min_per_bc))

m = np.asarray(minimizers)
mbc = np.asarray(min_per_bc)

plt.figure(1)

n, bins, p = plt.hist(m, bins = range(0, max_m+int(max_m/10), int(max_m/10)))
plt.xticks(bins)
plt.title("Minimizers")
plt.xlabel("Minimizer value")
plt.ylabel("Frequency")

plt.savefig("m.png")

plt.figure(2)
n, bins, p = plt.hist(mbc, bins = range(0, max_mbc+1, int(max_mbc/10)))
locs, labels = plt.xticks()
plt.xticks(bins)
plt.title("Minimizers per barcode")
plt.xlabel("Number of minimizers")
plt.ylabel("Frequency")

plt.savefig("mperb.png")
