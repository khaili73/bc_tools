#!/usr/bin/env python2.7

import matplotlib
matplotlib.use("AGG")

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle


#tiki
wm_file = "/home/kwegrzyn/data/comparison/all_white_mins.pickle"
plot_file = '/home/kwegrzyn/data/plots'

def undo_pickle(path):
    with open(path, 'rb') as f:
        obj = pickle.load(f)
    return obj


wm = undo_pickle(wm_file)
sns.distplot
plt.scatter(df.groups[:200], df.reads[:200])
plt.axis([0, 10, 0, 300])
plt.savefig('scatter_groups.pdf')
