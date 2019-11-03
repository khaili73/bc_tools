#!/usr/bin/env python3

import matplotlib
matplotlib.use("AGG")
import pandas as pd
import seaborn as sns

f = "/home/kwegrzyn/data/footprint_S_fixpvalue.tsv"

footprint_data = pd.read_csv(f, sep="\t", names = ["barcode", "minimizers_group"], name="minimizers number per barcode in S footprint group")
ax = sns.distplot(footprint_data, bins=20, kde=False)
sns_plot = ax.get_figure()
sns_plot.savefig("mins_number_per_bc_dist.pdf")
