#!/usr/bin/env python3
import matplotlib
matplotlib.use("AGG")
import seaborn as sns
import pickle


#tiki
# wm_file = "/home/kwegrzyn/data/comparison/all_white_mins.pickle"
wm_file = "/home/kwegrzyn/data/comparison/all_black_mins.pickle"
plot_file = '/home/kwegrzyn/data/plots'

def undo_pickle(path):
    with open(path, 'rb') as f:
        obj = pickle.load(f)
    return obj


wm = undo_pickle(wm_file)
ax = sns.distplot(wm, bins=20, kde=False)
sns_plot = ax.get_figure()
sns_plot.savefig("black_mins_dist.pdf")
