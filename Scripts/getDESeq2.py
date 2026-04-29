#!/usr/bin/env python3

# ------------------------------
# Script name: getDESeq2.py
# ------------------------------
# Purpose:
# Perform differential gene expression analysis using DESeq2 (via PyDESeq2)
# on RNA-seq count data. The script supports both single-factor designs
# (e.g., condition) and two-factor designs (e.g., genotype + treatment).

# How to run:
# python getDESeq2.py <counts.tsv> <metadata.csv> <output_dir> <prefix>

# Input:
# - counts.tsv: featureCounts output (cleaned), genes × samples
# - metadata.csv: sample metadata with experimental conditions
#	Examples for metadata file in folder "Metadata_examples"
# - output_dir: directory to store DESeq2 results
# - prefix: prefix for naming output files

# Output:
# - Differential expression results table(s) (log2FC, p-values, adjusted p-values)

# Notes:
# - Automatically detects experimental design from metadata columns
# - Applies appropriate contrasts depending on dataset conditions


import sys
import pandas as pd
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# ------------------------------
# Arguments
# ------------------------------
# Check that the script received the expected number of arguments
if len(sys.argv) != 5:
	sys.exit(
		"Usage: getDESeq2.py <counts.tsv> <metadata.csv> <outdir> <prefix>"
	)

# Read input paths from the command line
counts_file = sys.argv[1]
meta_file	= sys.argv[2]
outdir		= Path(sys.argv[3])
prefix		= sys.argv[4]

# Create the output directory if it does not already exist
outdir.mkdir(parents=True, exist_ok=True)

# ------------------------------
# Load counts (featureCounts format)
# ------------------------------
# Read the cleaned featureCounts table
counts = pd.read_csv(
	counts_file,
	sep="\t",
	comment="#"
)

# Set gene IDs as the row index
counts = counts.set_index("Geneid")

# Remove featureCounts metadata columns so only sample counts remain
counts = counts.drop(
	columns=["Chr", "Start", "End", "Strand", "Length"],
	errors="ignore"
)

# ------------------------------
# Load metadata
# ------------------------------
# Read the metadata file containing sample names and experimental labels
meta = pd.read_csv(
	meta_file,
	index_col=0
)

# Check that all metadata samples are present in the count matrix
missing = set(meta.index) - set(counts.columns)
if missing:
	raise ValueError(f"Samples in metadata not found in counts: {missing}")

# Reorder count columns to match metadata order, then transpose to samples × genes
counts = counts[meta.index].T

# ------------------------------
# Run DESeq2 with appropriate design
# ------------------------------
# Identify the metadata columns available for this dataset
meta_cols = set(meta.columns)

# ------------------------------
# Two-factor design: genotype + treatment
# ------------------------------
if {"genotype", "treatment"}.issubset(meta_cols):

	print("Detected 2×2 design")

	# Build the DESeq2 dataset using genotype and treatment plus interaction term
	dds = DeseqDataSet(
		counts=counts,
		metadata=meta,
		design="~ genotype + treatment + genotype:treatment"
	)

	# Run DESeq2 normalization and model fitting
	dds.deseq2()

	# Define the contrast of interest for this design
	stats = DeseqStats(
		dds,
		contrast=("treatment", "diABZI", "Vehicle")
	)

	# Generate the statistical summary for the contrast
	stats.summary()

	# Store the full results table
	res = stats.results_df

	# Save the DESeq2 results to a tsv file
	res.to_csv(
		outdir / f"{prefix}_DESeq2_WT_diABZI_vs_Vehicle.tsv",
		sep="\t"
	)

# ------------------------------
# Single-factor design: condition
# ------------------------------
elif {"condition"}.issubset(meta_cols):

	print("Detected single-factor design")

	# Build the DESeq2 dataset using condition as the only design variable
	dds = DeseqDataSet(
		counts=counts,
		metadata=meta,
		design="~ condition"
	)

	# Run DESeq2 normalization and model fitting
	dds.deseq2()

	# Collect the experimental conditions present in the metadata
	conditions = set(meta["condition"])
	print("Conditions detected:", conditions)

	# ------------------------------
	# Contrast logic
	# ------------------------------
	# Choose the contrast automatically based on the labels found in the metadata

	if {"STING WT", "STING KO"}.issubset(conditions):
		contrast = ("condition", "STING WT", "STING KO")

	elif {"STING-N153S", "PARENTAL"}.issubset(conditions):
		contrast = ("condition", "STING-N153S", "PARENTAL")

	elif {"Par", "sgSTING"}.issubset(conditions):
		contrast = ("condition", "Par", "sgSTING")

	elif {"STING", "GFP"}.issubset(conditions):
		contrast = ("condition", "STING", "GFP")

	elif {"agonist", "control"}.issubset(conditions):
		contrast = ("condition", "agonist", "control")

	elif {"Par", "Tax-R"}.issubset(conditions):
		contrast = ("condition", "Tax-R", "Par")

	elif {"siC", "siMYC"}.issubset(conditions):
		contrast = ("condition", "siMYC", "siC")

	else:
		raise ValueError(
			f"Unknown condition labels in metadata: {conditions}"
		)

	print("Using contrast:", contrast)

	# Run statistical testing for the selected contrast
	stats = DeseqStats(dds, contrast=contrast)

	# Generate the statistical summary
	stats.summary()

	# Store the full results table.
	res = stats.results_df

	# Save the DESeq2 results to a tsv file
	res.to_csv(
		outdir / f"{prefix}_DESeq2_results.tsv",
		sep="\t"
	)

# ------------------------------
# Unsupported metadata format
# ------------------------------
else:
	raise ValueError(
		"Metadata must contain either 'condition' or both "
		"'genotype' and 'treatment'."
	)

# Print a completion message when the analysis finishes successfully
print("DESeq2 analysis completed successfully.")