#!/usr/bin/env python2.7

import matplotlib
matplotlib.use("AGG")

import seaborn as sns
import pickle


#tiki
wm_file = "/home/kwegrzyn/data/comparison/all_white_mins.pickle"
plot_file = '/home/kwegrzyn/data/plots'

def undo_pickle(path):
    with open(path, 'rb') as f:
        obj = pickle.load(f)
    return obj


wm = undo_pickle(wm_file)
plot = sns.distplot(wm, bins=20, kde=False)
plot.savefig("white_mins_dist.pdf")
