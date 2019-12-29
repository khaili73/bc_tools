#!/usr/bin/env python3
import matplotlib
matplotlib.use("AGG")
import seaborn as sns
import pickle
import pandas as pd

#tiki
data = "/home/kwegrzyn/data/comparison_ML_fixpvalue.tsv"
plot_file = '/home/kwegrzyn/data/'

df = pd.read_csv(data, sep="\t", names = ["pvalue"], usecols=[2])
pvalue_dist = sns.distplot(df.pvalue, bins=range(0,0.000003,0.0000005), kde=False, )
#pvalue_dist.set(xlim=(0,0.000002))
fig = pvalue_dist.get_figure()
fig.savefig("ML_pvalue_dist_small.pdf")
