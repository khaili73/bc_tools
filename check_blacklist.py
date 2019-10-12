#!/usr/bin/env python3

import pickle

black_list_path = '/home/kwegrzyn/data/comparison/black_list.pickle'
minimizers_path = '/home/kwegrzyn/data/mod24_pickles/all_mins.pickle'


bl_stats_path = '/home/kwegrzyn/data/black_list_stats.txt'
bm_path = '/home/kwegrzyn/data/all_black_mins.pickle'
wm_path = '/home/kwegrzyn/data/all_white_mins.pickle'


def undo_pickle(path):
    with open(path, 'rb') as f:
        obj = pickle.load(f)
    return obj


def do_pickle(path, obj):
    with open(path, "wb") as f:
        pickle.dump(obj, f)


mins = undo_pickle(minimizers_path)
bl = undo_pickle(black_list_path)

all_white_mins = []
all_black_mins = []
for m in mins:
    if m in bl:
        all_black_mins.append(m)
    else:
        all_white_mins.append(m)


do_pickle(bm_path, all_black_mins)
do_pickle(wm_path, all_white_mins)


with open(bl_stats_path, 'w') as of:
    bl_len = str(len(bl))
    wl_len = str(len(set(all_white_mins)))
    bm_number = str(len(all_black_mins))
    wm_number = str(len(all_white_mins))
    of.write(f'white minimizers values number: {wl_len}')
    of.write(f'black minimizers values number: {bl_len}')
    of.write(f'all white mins number: {wm_number}')
    of.write(f'all black mins numebr: {bm_number}')
