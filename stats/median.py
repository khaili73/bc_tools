#!/usr/bin/env python3.6
import statistics

numbers = []
with open("/home/kwegrzyn/data/mins_per_bc.tsv", "r") as f:
    lines= f.readlines()
    for line in lines:
        n = line.strip().split("\t")[-1]
        numbers.append(int(n))
med = statistics.median(numbers)
print(med)
