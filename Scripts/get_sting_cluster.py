#!/usr/bin/env python3

# ------------------------------
# Script: get_sting_cluster.py
# ------------------------------
# Purpose:
# Identify a contiguous block of genes that follows the same expression pattern
# as STING in a heatmap based on z-scored log2TPM values

# The script:
# - reads a log2TPM matrix
# - filters it based on provided gene list
# - computes per-gene z-scores across samples
# - generates a clustered heatmap
# - extracts the largest contiguous block of genes with a consistent pattern
# - saves that gene block as a text file

# How to run:
# python get_sting_cluster.py <log2TPM.tsv> <gene_list.txt> <metadata.csv> <out_matrix.tsv> <out_heatmap.png>

# Input:
# - log2TPM.tsv: log2-transformed TPM matrix with a gene_name column
# - gene_list.txt: list of genes to test (1 gene per row) - this is the signature
# - metadata.csv: sample metadata with sample and condition columns
# - out_matrix.tsv: output file for the z-score matrix
# - out_heatmap.png: output file for the heatmap

# Output:
# - A text file containing the genes in the selected STING-like cluster


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from scipy.cluster.hierarchy import fcluster

# -----------------------------
# Argument checking
# -----------------------------
# Make sure the script received the expected number of command-line arguments
if len(sys.argv) != 6:
	print("Usage:")
	print("python get_sting_cluster.py <log2TPM.tsv> <gene_list.txt> <metadata.csv> <out_matrix.tsv> <out_heatmap.png>")
	sys.exit(1)

# Read input and output file paths from the command line
log2_file = sys.argv[1]
gene_list_file = sys.argv[2]
meta_file = sys.argv[3]
out_matrix = sys.argv[4]
out_heatmap = sys.argv[5]

# -----------------------------
# Load log2TPM
# -----------------------------
# Read the log2TPM matrix into a DataFrame
df = pd.read_csv(log2_file, sep="\t")

# Make sure the log2TPM file contains the gene_name column
if "gene_name" not in df.columns:
	raise ValueError("Column 'gene_name' not found in log2TPM file.")

# Use gene_name as the row index so genes can be selected by name
df = df.set_index("gene_name")

# -----------------------------
# Load gene list
# -----------------------------
# Read the input gene list into a Python list
gene_list = pd.read_csv(gene_list_file, header=None)[0].tolist()

# -----------------------------
# Load metadata
# -----------------------------
# Read sample metadata, which should contain sample names and conditions
meta = pd.read_csv(meta_file)

# Check that the metadata contains the required columns
if "sample" not in meta.columns or "condition" not in meta.columns:
	raise ValueError("Metadata must contain 'sample' and 'condition' columns.")

# Convert sample and condition columns to strings for consistency
meta["sample"] = meta["sample"].astype(str)
meta["condition"] = meta["condition"].astype(str)

# Extract the ordered list of sample names
samples = meta["sample"].tolist()

# -----------------------------
# Validate samples exist
# -----------------------------
# Confirm that all metadata samples are present in the log2TPM matrix
missing_samples = set(samples) - set(df.columns)
if missing_samples:
	raise ValueError(f"Samples not found in log2TPM file: {missing_samples}")

# -----------------------------
# Filter genes
# -----------------------------
# Keep only genes that are present in both the gene list and the expression matrix
present_genes = [g for g in gene_list if g in df.index]
missing_genes = [g for g in gene_list if g not in df.index]

# Save missing genes if needed for debugging or tracking
# missing_out = out_matrix.replace(".tsv", "_missing_genes.txt")
# with open(missing_out, "w") as f:
#	  for g in missing_genes:
#		  f.write(g + "\n")
# print(f"{len(missing_genes)} genes not found. Saved to {missing_out}")

# Stop the script if none of the genes in the list are present in the matrix
if len(present_genes) == 0:
	raise ValueError("None of the genes in the gene list were found.")

# -----------------------------
# Subset matrix
# -----------------------------
# Keep only the selected genes and the samples listed in the metadata
sub = df.loc[present_genes, samples]

# -----------------------------
# Compute Z-score (per gene)
# -----------------------------
# Compute the mean expression for each gene across samples
means = sub.mean(axis=1)

# Compute the standard deviation for each gene
# Replace zero standard deviation with NaN to avoid division errors
stds = sub.std(axis=1).replace(0, np.nan)

# Convert expression values to z-scores per gene
z = sub.sub(means, axis=0).div(stds, axis=0)

# -----------------------------
# Remove zero variance genes
# -----------------------------
# Identify genes that became all-NaN after z-scoring, usually because variance was zero
zero_var_genes = z.index[z.isna().all(axis=1)].tolist()

# Save zero-variance genes if needed for debugging
# zero_var_out = out_matrix.replace(".tsv", "_zero_variance_genes.txt")
# with open(zero_var_out, "w") as f:
#	  for g in zero_var_genes:
#		  f.write(g + "\n")
# print(f"{len(zero_var_genes)} zero-variance genes saved to {zero_var_out}")

# Remove genes with no variance from the clustering matrix
z_cluster = z.drop(zero_var_genes)

# -----------------------------
# Save Z-score matrix
# -----------------------------
# Create a row showing the sample conditions in the same sample order
condition_row = pd.DataFrame(
	[meta.set_index("sample").loc[samples, "condition"].values],
	columns=samples,
	index=["Gene"]
)

# Save the z-score matrix with the condition row if needed
# z_with_condition = pd.concat([condition_row, z])
# z_with_condition.to_csv(out_matrix, sep="\t")
# print(f"Z-score matrix saved to {out_matrix}")

# -----------------------------
# Generate heatmap
# -----------------------------
# Map each sample to its experimental condition for display
sample_to_condition = dict(zip(meta["sample"], meta["condition"]))

# Rename heatmap columns from sample IDs to condition labels
z_cluster.columns = [sample_to_condition[s] for s in samples]

# Set general plot appearance
sns.set(font_scale=0.7)

# Build a clustered heatmap using the z-score matrix
g = sns.clustermap(
	z_cluster,
	cmap="vlag",
	center=0,
	metric="euclidean",
	method="average",
	row_cluster=True,
	col_cluster=False,
	figsize=(8,10),
	cbar_pos=None
)

# Move the x-axis labels to the top of the heatmap
g.ax_heatmap.xaxis.tick_top()
g.ax_heatmap.xaxis.set_label_position('top')

# Save the heatmap if needed
# plt.savefig(out_heatmap, dpi=300)
plt.close()
# print(f"Heatmap saved to {out_heatmap}")

# -----------------------------
# Get dendrogram order
# -----------------------------
# Get the row order from the heatmap clustering dendrogram
row_order = g.dendrogram_row.reordered_ind

# Reorder genes according to the clustering result
ordered_genes = z_cluster.index[row_order]
ordered_matrix = z_cluster.loc[ordered_genes]

# -----------------------------
# Split samples into two groups
# -----------------------------
# Split the ordered matrix into two halves to compare expression patterns
n = ordered_matrix.shape[1]
group1 = ordered_matrix.iloc[:, :n//2]
group2 = ordered_matrix.iloc[:, n//2:]

# -----------------------------
# Compute difference score
# -----------------------------
# Compute the average expression difference between the two groups
gene_score = group2.mean(axis=1) - group1.mean(axis=1)

# Use the absolute value so the direction does not matter
pattern_mask = gene_score.abs() > 0.5

# -----------------------------
# Find largest contiguous block
# -----------------------------
# Track the longest consecutive region of genes that pass the pattern threshold
best_start = None
best_len = 0

current_start = None
current_len = 0

# Scan through the mask and identify the largest contiguous block
for i, val in enumerate(pattern_mask):

	if val:
		if current_start is None:
			current_start = i
			current_len = 1
		else:
			current_len += 1
	else:
		if current_len > best_len:
			best_len = current_len
			best_start = current_start
		current_start = None
		current_len = 0

# Check the final block at the end of the loop
if current_len > best_len:
	best_len = current_len
	best_start = current_start

# Extract the genes in the selected cluster block
block_genes = ordered_genes[best_start:best_start+best_len]

# -----------------------------
# Save gene list
# -----------------------------
# Save the selected genes to a text file, one gene per line
out_file = out_matrix.replace(".tsv", "_sting_cluster.txt")
pd.Series(block_genes).to_csv(out_file, index=False, header=False)

# Print a summary of the final result
print(f"{len(block_genes)} genes saved to {out_file}")
print()