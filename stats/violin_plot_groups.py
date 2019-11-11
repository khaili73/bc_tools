#!/usr/bin/env python3

import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

#local
#f = "./groups_stats.tsv"
#tiki
f = "/home/kwegrzyn/data/groups/groups_stats.tsv"

barcodes = pd.read_csv(f, sep="\t", usecols=[0,1,2,3], names = ["barcode", "read_number", "groups", "valid_read_number"])
#data.plot.scatter(x='groups', y='all_reads', c='DarkBlue')
#f, (ax) = plt.subplots(1, 1, figsize=(12, 4))
#f.suptitle('Clustering of barcodes reads mapping', fontsize=14)
#plt.xlabel('Groups')
#plt.ylabel('Number of reads')

violin_plot = sns.violinplot(x="groups", y="read_number",
               data=barcodes, inner="quart", linewidth=1.3)
#ax.set_xlabel("Number of groups",size = 12,alpha=0.8)
#ax.set_ylabel("Number of reads",size = 12,alpha=0.8)

violin_plot.set(ylim=(0,250),xlim=(0,6))
fig = violin_plot.get_figure()
fig.savefig('/home/kwegrzyn/data/groups/num_of_reads_per_group_violin_small.pdf')

