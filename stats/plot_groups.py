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
plt.scatter(barcodes['groups'], barcodes['read_number'],
            alpha=0.4, edgecolors='w')

plt.xlabel('Groups')
plt.ylabel('Number of reads')
plt.title('Clustering of barcodes reads mapping',y=1.05)


# Joint Plot
jp = sns.jointplot(x='groups', y='read_number', data=barcodes,
                   kind='reg', space=0, size=5, ratio=4)

plt.savefig('/home/kwegrzyn/data/groups/sns_scatter_groups.pdf')
