#!/usr/bin/env python3

import os
import pysam
from collections import defaultdict
import re

#locally
bam_dir = "/Users/khaili/mgr/test_data/sorted_bams/"
#tiki
#bam_dir = ""

with open("groups_stats.tsv", 'a') as f:
    #groups = defaultdict(list)
    for filename in os.listdir(bam_dir):
        groups = []
        chromosome_groups = defaultdict(list)
        samfile = pysam.AlignmentFile(bam_dir+filename, 'rb')
        bc = re.search("BX:Z:[ATCG]{16}-1", filename).group()
        reads_number = 0
        for read in samfile:
            reads_number += 1
            r_dict = read.to_dict()
            ch = r_dict['ref_name']
            chromosome_groups[ch].append(read)
        groups_number = 0
        valid_reads_number = 0
        for c in chromosome_groups:
            if len(chromosome_groups[c]) >= 10:
                subgroup = []
                last_position = chromosome_groups[c][0].reference_start
                for read in chromosome_groups[c]:
                    position = read.reference_start
                    if abs(position - last_position) < 10000:
                        subgroup.append(read)
                    last_position = read.reference_end
                if len(subgroup) >= 10:
                    #groups[bc].append(subgroup)
                    groups.append(subgroup)
        #groups_number = len(groups[bc])
        groups_number = len(groups)
        #valid_reads_number = sum([len(g) for g in groups[bc]])
        valid_reads_number = sum([len(g) for g in groups])
        f.write('\t'.join([str(x) for x in (bc, reads_number, groups_number, valid_reads_number, '\n')]))
#print(groups)
