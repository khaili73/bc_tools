#!/bin/bash

for file in /home/kwegrzyn/data/grouped/*
{
    bc=$(echo $file | tr "/" "." | tr "." "\n")
    set -- $bc
    bowtie2 --no-unal -p 4 -x genome -U $file -S "/home/kwegrzyn/data/index/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/sams/$5.sam" &
    while  [ `ps -ef | grep bowtie2 | grep -v grep |wc -l` -gt 10 ]
    do
        sleep 1;
    done
}
