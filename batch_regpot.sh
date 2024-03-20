#!/bin/bash

# Remove the file if it exists
if [ -f rp_score_all_TF.csv ]; then
    rm rp_score_all_TF.csv
fi

while IFS=, read -r gene_symbol species half_decay
do
    ./regpot_single_gene_hg38_mm10.py "$gene_symbol" "$species" "$half_decay"
    done < genes_oi.txt
