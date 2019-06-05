#!/usr/bin/env python3

import matplotlib
matplotlib.use("AGG")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#local
#f = "./groups_stats.tsv"
#tiki
f = "/home/kwegrzyn/data/groups_stats.tsv"

df = pd.read_csv(f, sep="\t", names = ["valid_reads", "groups", "reads", "NN"])
#data.plot.scatter(x='groups', y='all_reads', c='DarkBlue')
plt.scatter(df.groups, df.reads)
plt.savefig('scatter_groups.pdf')
