import pandas as pd
import numpy as np
import sys
import os

# ------------------------------
# Script: twoTPMmatrix.py
# ------------------------------
# Purpose:
# Convert raw gene counts into TPM (Transcripts Per Million) and log2(TPM + 1)
# to normalize for both sequencing depth and gene length

# How to run:
# python twoTPMmatrix.py <input_counts.tsv> <output_TPM.tsv> <output_log2TPM.tsv>

# Example:
# python3 twoTPM_matrix.py \
#	/N/slate/cma16/counts/GSE165825/all_counts_GSE165825_cleaned.txt \
#	/N/slate/cma16/counts/GSE165825/all_TPM_GSE165825.tsv \
#	/N/slate/cma16/counts/GSE165825/all_log2TPM_GSE165825.tsv

# Input:
# - Tab-separated counts file (featureCounts output, cleaned)
# - Must include "Geneid" and "Length" columns

# Output:
# - TPM matrix
# - log2(TPM + 1) matrix (used for downstream analysis and visualization)

# ------------------------------
# Argument Checking
# ------------------------------
if len(sys.argv) != 4:
	print("Usage: python twoTPM_matrix.py <input_counts> <output_TPM> <output_log2TPM>")
	sys.exit(1)

infile = sys.argv[1]
tpm_outfile = sys.argv[2]
log2_outfile = sys.argv[3]

if not os.path.exists(infile):
	print(f"ERROR: Input file '{infile}' does not exist.")
	sys.exit(1)

print(f"Reading counts from: {infile}")

# ------------------------------
# Read counts file
# ------------------------------
df = pd.read_csv(infile, sep="\t", comment="#")

# Required metadata columns from featureCounts
required_cols = ['Geneid', 'Length']
for col in required_cols:
	if col not in df.columns:
		print(f"ERROR: Required column '{col}' not found in input file.")
		sys.exit(1)

# Identify metadata and sample columns
meta_cols = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']
sample_cols = [col for col in df.columns if col not in meta_cols]

if len(sample_cols) == 0:
	print("ERROR: No sample columns detected in input file.")
	sys.exit(1)

# ------------------------------
# Clean sample names
# ------------------------------
clean_names = [
	os.path.basename(col).replace("_Aligned.sortedByCoord.out.bam", "")
	for col in sample_cols
]

rename_map = dict(zip(sample_cols, clean_names))
df.rename(columns=rename_map, inplace=True)

# ------------------------------
# Prepare data
# ------------------------------
length_kb = df['Length'] / 1000

# Avoid division by zero
length_kb.replace(0, np.nan, inplace=True)

counts = df[clean_names]

# Ensure counts are numeric
counts = counts.apply(pd.to_numeric, errors='coerce')

# ------------------------------
# TPM Calculation
# ------------------------------
# Step 1: Reads Per Kilobase (RPK)
rpk = counts.div(length_kb, axis=0)

# Step 2: Scaling factor (per sample)
scaling_factors = rpk.sum(axis=0) / 1e6

if (scaling_factors == 0).any():
	print("WARNING: One or more samples have zero scaling factor (possible empty counts).")

# Step 3: TPM normalization
tpm = rpk.div(scaling_factors, axis=1)

# Add Gene IDs back
tpm.insert(0, 'Geneid', df['Geneid'])

# ------------------------------
# Save TPM
# ------------------------------
tpm.to_csv(tpm_outfile, sep="\t", index=False)
print(f"TPM matrix saved to: {tpm_outfile}")

# ------------------------------
# log2(TPM + 1)
# ------------------------------
# Log-transform for downstream analysis (e.g., heatmaps)
log2_tpm = np.log2(tpm.iloc[:, 1:] + 1)
log2_tpm.insert(0, 'Geneid', tpm['Geneid'])

log2_tpm.to_csv(log2_outfile, sep="\t", index=False)
print(f"log2(TPM+1) matrix saved to: {log2_outfile}")