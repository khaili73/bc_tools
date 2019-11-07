#!/usr/bin/env python3

import os
from Bio import SeqIO

path = ''
directory = os.fsencode(path)

for file in os.listdir(directory):
     filename = os.fsdecode(file)
     if filename.endswith(".fastq"):
         # print(os.path.join(directory, filename))
         for record in SeqIO.parse(os.path.join(directory, filename), "fastq"):
             print(record.id)
         continue
     else:
         continue
