#!/bin/bash

for file in /home/kwegrzyn/data/grouped/*
{
    bc=$(echo $file | tr "/" "." | tr "." "\n")
    set -- $bc
    bowtie2 --no-unal -x genome -S "$5.sam" &
    while  [ `ps -ef | grep bowtie2 | grep -v grep |wc -l` -gt 14 ]
    do
        sleep 1;
    done
}
