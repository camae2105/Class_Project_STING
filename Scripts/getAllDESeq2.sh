#!/bin/bash

#SBATCH --job-name=deseq_all_datasets
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cma16@iu.edu	# Add email here!
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=/N/slate/cma16/logs/%x_%j.out	# CHANGE user cma16
#SBATCH --error=/N/slate/cma16/logs/%x_%j.err	# CHANGE user cma16
#SBATCH -A r00235

# Exit if any command fails, undefined variable is used, or pipe fails
set -euo pipefail

# Load Python for running the DESeq2 wrapper script
module load python

# Define directories 
BASE=/N/slate/cma16 # Slate directory - CHANGE the user cma16 here!!!
COUNTS_BASE=$BASE/counts1 # Dir with count matrices for each dataset
METADATA_DIR=$BASE/metadata # Dir with metadata files for each dataset/experiment
OUT_BASE=$BASE/deseq # Output dir for DESeq2 results
SCRIPT_DIR=$BASE/python_scripts # Python scripts dir

# List datasets you want to process
# Specify the datasets to loop through in this batch job
DATASETS=(
	GSE165825
	GSE166209
	GSE196323
	GSE271679
	GSE283768
)

# Loop through each dataset and run the DESeq2 workflow
for DATASET in "${DATASETS[@]}"
do
	echo "=============================="
	echo "Processing dataset: $DATASET"
	echo "=============================="

	# Path to the cleaned counts file for this dataset
	COUNTS_FILE=$COUNTS_BASE/$DATASET/all_counts_${DATASET}_cleaned.txt

	# Output directory for this dataset
	OUTDIR=$OUT_BASE/$DATASET

	# Create the dataset-specific output directory if needed
	mkdir -p "$OUTDIR"

	# Skip the dataset if the counts file is missing
	if [ ! -f "$COUNTS_FILE" ]; then
		echo "WARNING: Counts file not found for $DATASET — skipping"
		continue
	fi

	# Ensure wildcard returns empty list if no files match (avoids fake filenames)
	shopt -s nullglob

	# Collect all metadata files for the current dataset
	# Metadata examples in folder "Metadata_examples"
	meta_files=($METADATA_DIR/metadata_${DATASET}_*.csv)

	# Skip the dataset if no metadata files are available
	if [ ${#meta_files[@]} -eq 0 ]; then
		echo "WARNING: No metadata files found for $DATASET — skipping"
		continue
	fi

	# Run DESeq2 once for each metadata file/experiment
	for meta in "${meta_files[@]}"
	do
		# Use the metadata filename as the output prefix
		prefix=$(basename "$meta" .csv)

		echo "Running DESeq2 for: $prefix"

		# Run the custom Python script that performs DESeq2 analysis
		python "$SCRIPT_DIR/getDESeq2.py" \
			"$COUNTS_FILE" \
			"$meta" \
			"$OUTDIR" \
			"$prefix"
	done
done

# Print a completion message when all datasets are finished
echo "All DESeq2 analyses completed."