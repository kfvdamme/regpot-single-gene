#!/usr/bin/env python

import argparse
import itertools
import sqlite3

from pybedtools import BedTool
import pandas as pd

# essential data locations
gene_database = 'regpot_single_gene.db'
hg38_remap2022 = 'remap2022/remap2022_nr_macs2_hg38_v1_0.bed'
mm10_remap2022 = 'remap2022/remap2022_nr_macs2_mm10_v1_0.bed'


def get_gene_info(gene_symbol, species, conn):
    
    c = conn.cursor()

    # Determine the table name
    table_name = species + '_ncbiRefSeq'
    # Query the database
    c.execute(f"SELECT chrom, strand, txStart, txEnd, name2 FROM {table_name} WHERE name2=?", (gene_symbol,))

    # Fetch the results
    results = c.fetchone()

    # Calculate TSS
    chrom, strand, txStart, txEnd, name2 = results
    TSS = txStart if strand == '+' else txEnd

    
    # Return the results
    return chrom, strand, txStart, txEnd, name2, TSS, species


def drop_table(table_name, conn):
    # Connect to the SQLite database
    cursor = conn.cursor()

    # Execute the DROP TABLE statement
    cursor.execute(f"DROP TABLE IF EXISTS {table_name}")

    # Commit the changes and close the connection
    conn.commit()
    #conn.close()

def create_and_insert_into_gene_info(gene_info, conn):
    # Create a new table
    c = conn.cursor()
    c.execute("""
        CREATE TABLE IF NOT EXISTS gene_info (
            chrom TEXT,
            strand TEXT,
            txStart INTEGER,
            txEnd INTEGER,
            name2 TEXT,
            TSS INTEGER,
            species TEXT
        )
    """)

    # Insert the gene info into the table
    c.execute("""
        INSERT INTO gene_info (chrom, strand, txStart, txEnd, name2, TSS, species)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    """, gene_info)

    # Commit the changes
    conn.commit()

def get_overlapping_regions(gene_name, half_decay, conn):
    # 1. Extract the TSS and species of the gene from the gene_info table
    c = conn.cursor()
    c.execute("SELECT chrom, strand, TSS, species FROM gene_info WHERE name2=?", (gene_name,))
    chrom, strand, TSS, species = c.fetchone()

    # 2. Create a region from the TSS by extending it to the left and right with the half_decay distance
    half_decay_distance = int(half_decay[:-2]) * 1000  # Convert from kb to base pairs
    region_start = TSS - half_decay_distance
    region_end = TSS + half_decay_distance

    # Create a region in the specified format
    region = {'chrom': chrom, 'strand': strand, 'start': region_start, 'end': region_end, 'TSS': TSS}

    # 3. Determine the bed file to use based on the species
    bed_file = hg38_remap2022 if species == 'hg38' else mm10_remap2022

    # Create a BedTool object for the region and the bed file
    region_bed = BedTool(f"{region['chrom']} {region['start']} {region['end']}", from_string=True)
    bed = BedTool(bed_file)

    # Get the overlapping regions
    overlapping_regions = bed.intersect(region_bed)

    return overlapping_regions

def create_table_from_overlapping_regions(overlapping_regions, species, conn):
    # Create a new table
    c = conn.cursor()
    c.execute("""
        CREATE TABLE IF NOT EXISTS overlapping_regions (
            chrom TEXT,
            start INTEGER,
            end INTEGER,
            name TEXT,
            peak INTEGER,
            TF TEXT,
            species TEXT
        )
    """)

    # Parse the overlapping_regions and insert the data into the table
    for region in overlapping_regions:

        chrom = region.chrom
        start = region.start
        end = region.end
        name = region.name
        peak = region[6]  # Access the thickStart field
        TF = name.split(':')[0]

        c.execute("""
            INSERT INTO overlapping_regions (chrom, start, end, name, peak, TF, species)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """, (chrom, start, end, name, peak, TF, species))

    # Commit the changes
    conn.commit()

def list_TFs(table_name, species, conn):
    # Connect to the SQLite database
    cursor = conn.cursor()

    # Execute the query
    cursor.execute(f"SELECT DISTINCT TF FROM {table_name} WHERE species = '{species}'")

    # Fetch the results
    results = cursor.fetchall()

      # Extract the TFs from the results and return them
    TFs = [row[0] for row in results]
    return TFs

def get_peak_TSS_distances(species, TF, conn):
    # Get a cursor
    c = conn.cursor()

    # Perform the JOIN operation and calculate the distances
    c.execute("""
        SELECT ABS(overlapping_regions.peak - gene_info.TSS) AS distance
        FROM overlapping_regions
        JOIN gene_info ON overlapping_regions.chrom = gene_info.chrom AND overlapping_regions.species = gene_info.species
        WHERE overlapping_regions.species = ? AND overlapping_regions.TF = ?
    """, (species, TF))

    # Fetch all the distances
    distances = [row[0] for row in c.fetchall()]

    return distances

def calculate_fraction(distances):
    # Check if distances is empty
    if not distances:
        return 0

    small_distances = [d for d in distances if d <= 1000]

    # Calculate the fraction
    fraction = len(small_distances) / len(distances)

    return fraction

def calculate_rp_score(d0, distances):
    rp_score = 0
    for distance in distances:
        if distance / d0 > 100:  # Adjust the threshold as needed
            weight = 0
        else:
            weight = 1 / (2 ** (distance / d0))
        rp_score += weight
    return rp_score

def calculate_regulatory_potential(TF, species, conn):
    # Get the distances
    distances = get_peak_TSS_distances(species, TF, conn)

    # Calculate the fraction
    fraction = calculate_fraction(distances)

    # Set d0 based on the fraction
    d0 = 1000 if fraction > 0.2 else 10000

    # Calculate the regulatory potential score
    rp_score = calculate_rp_score(d0, distances)

    return fraction, rp_score

def create_rp_score_dataframe(TFs, gene_symbol, species, conn):
    # Initialize a list to store the data
    data = []

    # Calculate the rp_score for each TF and species
    for TF in TFs:
        rp_score = calculate_regulatory_potential(TF, species, conn)[1] 
        #data.append([gene_symbol, TF, rp_score, species])
        data.append([TF, gene_symbol, species, rp_score])

    # Create a DataFrame from the data
    df = pd.DataFrame(data, columns=['TF', 'Gene', 'species', 'RP_score'])

    return df


def main(gene_symbol, species, half_decay):
    print(f"Gene symbol: {gene_symbol}")
    print(f"Species: {species}")
    print(f"Half-decay distance: {half_decay}")

  # Connect to the SQLite database
    conn = sqlite3.connect(gene_database)

    gene_info = get_gene_info(gene_symbol, species, conn)
    print(gene_info)
    
    # remove any existing gene_info table
    drop_table('gene_info', conn)
    
    # Create the gene_info table and insert the gene info
    create_and_insert_into_gene_info(gene_info, conn)


    overlapping_regions = get_overlapping_regions(gene_symbol, half_decay, conn)
    
    # remove any existing overlapping_regions table
    drop_table('overlapping_regions', conn)

    create_table_from_overlapping_regions(overlapping_regions, species, conn)

    TFs = list_TFs('overlapping_regions', species, conn)
    print(f"Transcription factors: {TFs}")

    rp_score_all_TF_df = create_rp_score_dataframe(TFs, gene_symbol, species, conn)
    #print(rp_score_all_TF_df)

    with open('rp_score_all_TF.csv', 'a') as f:
        rp_score_all_TF_df.to_csv(f, header=False, index=False)

    # Close the connection
    conn.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Process a single gene - human or mouse.",
                                     usage="python regpot_gene_hg38_mm10.py <gene_symbol> <species> <half_decay>")
    parser.add_argument("gene_symbol", help="The gene symbol")
    parser.add_argument("species", help="species of the gene, either 'hg38' or 'mm10'")
    parser.add_argument("half_decay", choices=["1kb", "10kb", "100kb"], help="The half-decay distance")

    args = parser.parse_args()

    main(args.gene_symbol, args.species, args.half_decay)
