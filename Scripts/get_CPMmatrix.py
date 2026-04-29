#!/usr/bin/env python3

# ------------------------------
# Script name: getCPMmatrix.py
# ------------------------------
# Purpose:
# Convert raw gene counts (featureCounts output) into CPM (Counts Per Million)
# to normalize for sequencing depth across samples

# How to run:
# python getCPMmatrix.py <input_counts.tsv> <output_cpm.tsv>

# Input:
# - Tab-separated counts file (featureCounts output, cleaned)
# - Must include Geneid column and sample count columns

# Output:
# - CPM-normalized matrix (same structure as input, normalized values)


# Import required libraries
# pandas → for handling tabular data
# sys → for reading command-line arguments
# os → for checking file existence
import pandas as pd
import sys
import os

# ------------------------------
# 1. Check command-line arguments
# ------------------------------
# The script expects exactly two arguments:
#	1) Input counts file
#	2) Output CPM file
if len(sys.argv) != 3:
	print("Usage: python compute_cpm.py <input_counts.tsv> <output_cpm.tsv>")
	sys.exit(1)

# Assign arguments to variables
infile = sys.argv[1]
outfile = sys.argv[2]

# Check that the input file exists
if not os.path.exists(infile):
	print(f"ERROR: Input file '{infile}' not found.")
	sys.exit(1)

# ------------------------------
# 2. Load the counts matrix
# ------------------------------
# We assume this is the cleaned featureCounts file:
#	- Rows = genes
#	- Columns = metadata + sample counts
#	- Tab-separated
df = pd.read_csv(infile, sep="\t", comment="#")

print(f"Loaded counts file: {infile}")
print(f"Shape of counts matrix: {df.shape}")

# ------------------------------
# 3. Identify sample columns
# ------------------------------
# featureCounts output usually includes these metadata columns:
meta_cols = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']

# Sample columns are everything not in metadata
sample_cols = [col for col in df.columns if col not in meta_cols]

if len(sample_cols) == 0:
	print("ERROR: No sample columns detected in input file.")
	sys.exit(1)

print(f"Detected {len(sample_cols)} samples.")

# ------------------------------
# 4. Extract raw counts only
# ------------------------------
# Create a new DataFrame containing only count values
counts = df[sample_cols]

# Convert to numeric in case there are formatting issues
counts = counts.apply(pd.to_numeric, errors='coerce')

# ------------------------------
# 5. Compute library sizes
# ------------------------------
# Library size = total reads per sample
# This is used to normalize counts
lib_sizes = counts.sum(axis=0)

# Check for samples with zero library size (should never happen)
if (lib_sizes == 0).any():
	print("WARNING: One or more samples have total count = 0.")

# ------------------------------
# 6. Compute CPM (Counts Per Million)
# ------------------------------
# CPM formula:
#	CPM = (raw_count / total_reads_in_sample) * 1,000,000
cpm = counts.div(lib_sizes, axis=1) * 1e6

# ------------------------------
# 7. Add Gene IDs back to matrix
# ------------------------------
# Insert Geneid column at first position
cpm.insert(0, "Geneid", df["Geneid"])

# ------------------------------
# 8. Save CPM matrix
# ------------------------------
cpm.to_csv(outfile, sep="\t", index=False)

print(f"CPM matrix saved to: {outfile}")