#!/usr/bin/env python2.7

import matplotlib
matplotlib.use("AGG")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Plot histogram for one column of data')
parser.add_argument('infiles', metavar=('input'), help='input tsv files', nargs='+')
parser.add_argument('--out', default=os.environ['HOME'], help='output files directory')
args = parser.parse_args()

for f in args.infiles:
    name = f.split(".")[0]
    print(name)
    data = pd.read_csv(f, sep="\t", names = ["Minimizers"])
    stats = data.describe().to_csv(args.out + "/%s_description.tsv" % name, sep="\t")
    #max_m = data.max()
    max_m = 1000
    data.plot.hist(bins= range(0, max_m+max_m/10, max_m/10))
    plt.title('Minimizers occurance')
    plt.xlabel('Minimizers values')
    plt.grid(axis='y', linestyle='--')
    plt.savefig(args.out + "/%s_hist1000.png" % name)
