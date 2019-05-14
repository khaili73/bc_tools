#!/usr/bin/env python3

import pysam
import os

dir_str = "/Users/khaili/mgr/test_data/sams/"
directory = os.fsencode(dir_str)

for f in os.listdir(directory):
    filename = os.fsdecode(f)
    if filename.endswith(".sam"):
        file_path = os.path.join(dir_str, filename)
        samfile = pysam.AlignmentFile(file_path, "r")

