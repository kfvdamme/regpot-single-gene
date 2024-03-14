This Python script calculates the regulatory potential (RP) score for a given gene in either the human (hg38) or mouse (mm10) genome. The RP score is a measure of the likelihood that a transcription factor (TF) regulates a gene, based on the proximity of the TF's binding sites to the gene's transcription start site (TSS).

Here's a step-by-step breakdown of what the script does:

1. It takes three command-line arguments: the gene symbol, the species (either 'hg38' or 'mm10'), and the half-decay distance (either '1kb', '10kb', or '100kb').

2. It connects to a SQLite database that contains information about the genes and their transcription start sites (TSS).

3. It retrieves the gene's information from the database, including its chromosomal location, strand, transcription start and end sites, and species.

4. It creates a new table in the database to store the gene's information.

5. It calculates a region around the gene's TSS, extending to the left and right by the half-decay distance.

6. It determines which bed file to use based on the species, and creates a BedTool object for the region and the bed file.

7. It finds the regions in the bed file that overlap with the gene's region.

8. It creates a new table in the database to store the overlapping regions.

9. It retrieves a list of distinct transcription factors (TFs) from the overlapping regions.

10. For each TF, it calculates the distances between the TF's peak binding sites and the gene's TSS.

11. It calculates the fraction of distances that are less than or equal to 1000 base pairs.

12. It sets a parameter `d0` based on the fraction, and calculates the RP score for each TF.

13. It creates a DataFrame that contains the RP scores for all TFs, and appends it to a CSV file.

14. Finally, it closes the connection to the database.
