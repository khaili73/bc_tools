#!/usr/bin/env python3

import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

#local
#f = "./groups_stats.tsv"
#tiki
f = "/home/kwegrzyn/data/groups_stats.tsv"

barcodes = pd.read_csv(f, sep="\t", names = ["barcode", "read_number", "groups", "valid_read_number"])
#data.plot.scatter(x='groups', y='all_reads', c='DarkBlue')
f, (ax) = plt.subplots(1, 1, figsize=(12, 4))
f.suptitle('Clustering of barcodes reads mapping', fontsize=14)
plt.xlabel('Groups')
plt.ylabel('Number of reads')

sns.violinplot(x="groups", y="read_number",
               data=barcodes, inner="quart", linewidth=1.3, ax=ax)
ax.set_xlabel("Number of groups",size = 12,alpha=0.8)
ax.set_ylabel("Number of reads",size = 12,alpha=0.8)

plt.savefig('/home/kwegrzyn/data/groups/sns_scatter_groups.pdf')
