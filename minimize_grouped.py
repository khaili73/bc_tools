#!/usr/bin/env python3
import os
import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import namedtuple

modulo = 24
minimizers_values_file = 'all_mins.tsv'
barcode_minimizers_file = 'bc_mins.tsv'
grouped_path = '/home/kwegrzyn/data/grouped/'


def minimize(k, sequence):
    dna_strand = Seq(sequence, generic_dna)
    rc_sequence = dna_strand.reverse_complement()
    minimizer = abs(hash(str(sequence)[:3]))
    for seq in str(sequence), str(rc_sequence):
        for i in range(len(seq) - k + 1):
            hashed_kmer = hash(seq[i:i + k])
            if modulo:
                hashed_kmer = hashed_kmer % modulo  # negative values
            if abs(hashed_kmer) <= minimizer:  # choose rigthmost lowest value
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


with open(minimizers_values_file) as min_file, open(barcode_minimizers_file) as bc_file:
    directory = os.fsencode(grouped_path)
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".fastq"):
            with open(filename, 'r') as fastq:
                records_parser = parse_fastq(fastq)
                mins = []
                for record in records_parser:
                    min = minimize(record.sequence)
                    mins.append(min)
                    min_file.write(min)
                bc = get_barcode(filename)
                bc_file.write(bc, len(mins), mins)
