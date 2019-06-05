#!/bin/bash

#for local testing
#sams=/Users/khaili/mgr/test_data/sams/
#bams=/Users/khaili/mgr/test_data/bams/
#sorted=/Users/khaili/mgr/test_data/sorted_bams/
sams=/home/kwegrzyn/data/sams/
bams=/home/kwegrzyn/data/bams/
sorted=/home/kwegrzyn/data/sorted_bams/

for file in $sams*
{
    bc=$(echo $file | tr "/" "." | tr "." "\n")
    set -- $bc
    samtools view -S -b $sams$5.sam > $bams$5.bam
    samtools view -H $bams$5.bam > header.sam
    #filter out unmapped [-F4] and mapped multiple times [XS:] reads
    samtools view -F4 $bams$5.bam | grep -v "XS:" | cat header.sam - | samtools view -b | samtools sort - -o $sorted$5.bam
    rm header.sam
}
