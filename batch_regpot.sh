#!/bin/bash

while IFS=, read -r gene_symbol species
do
    ./regpot_single_gene_hg38_mm10.py "$gene_symbol" "$species" 10kb
done < genes_oi.txt
