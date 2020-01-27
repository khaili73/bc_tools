#!/usr/bin/env python3
import os
import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import namedtuple

modulo = 2**24
minimizers_values_file = 'all_mins.tsv'
barcode_minimizers_file = 'bc_mins.tsv'
#grouped_path = '/home/kwegrzyn/data/grouped/'
grouped_path = '/Users/khaili/mgr/bc_tools/grouped/'


def minimize(k, sequence):
    dna_strand = Seq(sequence, generic_dna)
    rc_sequence = dna_strand.reverse_complement()
    minimizer = abs(hash(str(sequence)[:k])) % modulo
    print("initial: ", minimizer)
    for seq in str(sequence), str(rc_sequence):
        for i in range(len(seq) - k + 1):
            print(seq[i:i + k])
            hashed_kmer = hash(seq[i:i + k]) % modulo #NEGATIVE VALUES
            print("hashed: ", hashed_kmer)
            if abs(hashed_kmer) < minimizer:  # choose rigthmost lowest value
                minimizer = abs(hashed_kmer)
    return minimizer


def parse_fastq(handle):
    record = namedtuple('record', 'header, sequence, plusline, quality')
    while True:
        try:
            record_header = next(handle).strip("\n")
            record_seq = next(handle).strip("\n")
            record_plus = next(handle).strip("\n")
            record_quality = next(handle).strip("\n")
            yield record(record_header, record_seq, record_plus, record_quality)
        except StopIteration:
            handle.close()
            break


def get_barcode(record):
    if re.search("BX:Z:[ATCG]{16}-1", record.header):
        return re.search("BX:Z:[ATCG]{16}-1", record.header).group()
    else:
        return "BX:Z:0"


with open(minimizers_values_file, 'w') as min_file, open(barcode_minimizers_file, 'w') as bc_file:
    directory = os.fsencode(grouped_path)
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".fastq"):
            with open(grouped_path+filename, 'r') as fastq:
                records_parser = parse_fastq(fastq)
                mins = []
                for record in records_parser:
                    min = minimize(21, record.sequence)
                    mins.append(str(min))
                    min_file.write(str(min)+"\n")
                    bc = get_barcode(record)
                bc_file.write("\t".join((bc, str(len(mins)), str(",".join(mins))+"\n")))
