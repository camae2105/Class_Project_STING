import pandas as pd
import sys

# ------------------------------
# Script name: add_gene_names.py
# ------------------------------
# Purpose:
# Add gene names to a TPM matrix by merging a GeneID-to-name mapping file

# How to run:
# python add_gene_names.py <input_TPM.tsv> <gene_map.tsv> <output.tsv>


# Inputs
# 1) TPM matrix file (must contain "Geneid" column)
tpm_file = sys.argv[1]

# 2) Gene ID → gene name mapping file
gene_map_file = sys.argv[2]

# 3) Output file with gene names added
output_file = sys.argv[3]

# ------------------------------
# Load TPM matrix
# ------------------------------
df = pd.read_csv(tpm_file, sep="\t")

# ------------------------------
# Load gene mapping file
# ------------------------------
# Expecting two columns: Geneid and gene_name
gene_map = pd.read_csv(gene_map_file, sep="\t", header=None, names=["Geneid", "gene_name"])

# ------------------------------
# Merge TPM with gene names
# ------------------------------
# Left join ensures all genes in TPM are kept
df = gene_map.merge(df, on="Geneid", how="left")

# ------------------------------
# Save final output
# ------------------------------
df.to_csv(output_file, sep="\t", index=False)

print(f"Gene names added!!! Output saved to: {output_file}")