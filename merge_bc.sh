#!/bin/bash

awk -F'\t' -v OFS='\t' '{
    if ($1 in a) {
        a[$1] += $2;
    } else {
        a[$1] = $2;
    }
}
END { for (i in a) print i, a[i]}' < mins_per_bc_mod216.tsv
