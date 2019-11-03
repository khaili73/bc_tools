#!/usr/bin/env python3

import matplotlib
matplotlib.use("AGG")
import pandas as pd
import seaborn as sns

f = "/home/kwegrzyn/data/footprint_S_fixpvalue.tsv"

footprint_data = pd.read_csv(f, sep="\t", usecols=[0,1,2], names = ["barcode", "footprint", "minimizers"])
footprint_data['minimizers'] = footprint_data.minimizers.map(lambda x: [i for i in str(x).split(",")])
footprint_data['minimizers_len'] = footprint_data.minimizers.apply(len)
ax = sns.distplot(footprint_data['minimizers_len'], bins=range(1,31), kde=False)
sns_plot = ax.get_figure()
sns_plot.savefig("S_mins_number_per_bc_dist_30bins.pdf")
