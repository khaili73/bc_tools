#!/bin/bash

awk -F'\t' -v OFS='\t' '{
    if ($1 in a) {
        a[$1] += $3;
    } else {
        a[$1] = $3;
    }
}
END { for (i in a) print i, a[i]}' < "test.tsv"
