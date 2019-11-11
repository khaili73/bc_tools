#!/usr/bin/env python3

import gzip
import re
from collections import namedtuple

barcodes_dir = "/home/kwegrzyn/data/post-longranger"
reads = "/outs/barcoded.fastq.gz"
paths = [barcodes_dir + "/barcodes" + str(x) + reads for x in range(1,9)]

all_bc = []

def parse_fastq(path):
    record = namedtuple('record', 'header, sequence, quality')
    handle = open(path)
    while True:
        try:
            record_header = next(handle).strip("\n")
            record_seq = next(handle).strip("\n")
            next(handle) # plus line
            record_quality = next(handle).strip("\n")
            yield record(record_header, record_seq, record_quality)
        except StopIteration:
            handle.close()
            break

for path in paths:
    with gzip.open(path, 'rt') as file:
        read_parser = parse_fastq(file)
        for record in read_parser:
            if re.search("BX:Z:[ATCG]{16}-1", record.header):
                bc = re.search("BX:Z:[ATCG]{16}-1", record.header).group()
                all_bc.append(bc)

uniq_bc = set(all_bc)

with open("barcodes_count.txt", 'w') as ofile:
    ofile.write("Number of barcodes: {bc_num}\n".format(bc_num=len(uniq_bc)))
    ofile.write("Number of reads with barcode: {reads_num}\n".format(reads_num=len(all_bc)))