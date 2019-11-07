#!/usr/bin/env python3

import gzip
import re
from Bio import SeqIO

barcodes_dir = "/home/kwegrzyn/data/post-longranger"
reads = "/outs/barcoded.fastq.gz"
paths = [barcodes_dir + "/barcodes" + str(x) + reads for x in range(1,9)]

all_bc = []

for path in paths:
    with gzip.open(path, 'rt') as file:
        read_parser = SeqIO.parse(file, "fastq")
        for record in read_parser:
            if re.search("BX:Z:[ATCG]{16}-1", record.description):
                bc = re.search("BX:Z:[ATCG]{16}-1", record.description).group()
                all_bc.append(bc)

uniq_bc = set(all_bc)

with open("barcodes_count.txt", 'w') as ofile:
    ofile.write("Number of barcodes: {bc_num}\n".format(bc_num=len(uniq_bc)))
    ofile.write("Number of reads with barcode: {reads_num}\n".format(reads_num=len(all_bc)))