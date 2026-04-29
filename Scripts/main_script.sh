#!/bin/bash

# ---------------------------------------------------
# Main workflow script for STING signature project
# ---------------------------------------------------
# This script shows the order of the RNA-seq processing pipeline.
# Each step corresponds to a dataset-specific job script.

# Workflow:
# 1) Download SRA reads and convert them to FASTQ
# 2) Trim FASTQ files
# 3) Align reads with STAR
# 4) Generate counts, CPM, TPM, and log2TPM
# 5) Run differential expression analysis (DESeq2)

# The same pipeline was applied to all datasets.
# For larger datasets (GSE196323 and GSE271679),
# the alignment step was split into multiple jobs.
# ---------------------------------------------------

# ---------------------------------------------------
# Environment setup (run once before pipeline)
# ---------------------------------------------------

# Load Python module
module load python

# Install PyDESeq2 (Python implementation of DESeq2)
pip install --user pydeseq2
# Differential expression analysis was performed using PyDESeq2

# ---------------------------------------------------
# Software and dependencies
# ---------------------------------------------------
# The following tools were used in this pipeline:
#
# - SRA Toolkit (prefetch, fasterq-dump)
# - fastp (read trimming)
# - STAR (read alignment)
# - Subread (featureCounts for gene counting)
# - Python (pandas, numpy, seaborn, matplotlib, scipy, pydeseq2)
#
# Note:
# All tools except PyDESeq2 were available through HPC environment modules
# and loaded within each job script.
# ---------------------------------------------------

DATASETS=(
	GSE165825
	GSE166209
	GSE196323
	GSE271679
	GSE283768
)

for DATASET in "${DATASETS[@]}"
do
	echo "======================================"
	echo "Processing dataset: $DATASET"
	echo "======================================"

	# ---------------------------------------------------
	# STEP 1: Download SRA reads and convert to FASTQ
	# ---------------------------------------------------
	# Raw sequencing runs were downloaded from SRA and
	# converted into FASTQ format.
	sbatch get_fastq_${DATASET}.sh

	# ---------------------------------------------------
	# STEP 2: Trim FASTQ files
	# ---------------------------------------------------
	# FASTQ files were quality-trimmed to remove adapters
	# and low-quality bases.
	sbatch trim_fastq_${DATASET}.sh

	# ---------------------------------------------------
	# STEP 3: Align reads with STAR
	# ---------------------------------------------------
	# Trimmed reads were aligned to the human reference genome.
	
	# Most datasets use a single alignment job.
	# Larger datasets were split into multiple jobs to reduce runtime.

	if [[ "$DATASET" == "GSE196323" ]]; then
		# Alignment split into 3 jobs
		sbatch star_align1_${DATASET}.sh
		sbatch star_align2_${DATASET}.sh
		sbatch star_align3_${DATASET}.sh

	elif [[ "$DATASET" == "GSE271679" ]]; then
		# Alignment split into 2 jobs
		sbatch star_align1_${DATASET}.sh
		sbatch star_align2_${DATASET}.sh

	else
		# Standard single alignment job
		sbatch star_align_${DATASET}.sh
	fi

	# ---------------------------------------------------
	# STEP 4: Generate counts and normalized matrices
	# ---------------------------------------------------
	# Aligned reads were summarized to gene-level counts
	# using featureCounts. Counts were then normalized to:
	# - CPM (Counts Per Million)
	# - TPM (Transcripts Per Million)
	# - log2(TPM + 1) for downstream analysis
	sbatch getCountsTPM_${DATASET}.sh

	# ---------------------------------------------------
	# STEP 5: Differential expression analysis
	# ---------------------------------------------------
	# Cleaned count matrices and metadata were used to run
	# DESeq2 for each experiment within each dataset.
	sbatch getAllDESeq2.sh

	echo "Submitted all jobs for $DATASET"
	echo
done

echo "Workflow submission complete."

### Additional steps for differentially expressed genes analysis

# ---------------------------------------------------
# STEP 6: Extract STING-associated gene lists
# ---------------------------------------------------
# Differentially expressed genes were filtered by log2 fold change and
# p-value thresholds to generate STING-associated gene lists.

# Example:
# python getGeneList.py -i <DESeq2_results.tsv> --lfc 0 --pval 0.05

# ---------------------------------------------------
# STEP 7: Identify STING-like clusters in heatmaps
# ---------------------------------------------------
# Heatmap-based clustering was used to identify blocks of genes that
# followed the same expression pattern as STING.

# Example:
# python get_sting_cluster.py <log2TPM.tsv> <gene_list.txt> <metadata.csv> <out_matrix.tsv> <out_heatmap.png>
