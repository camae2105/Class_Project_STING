#!/bin/bash

############################
# SLURM JOB CONFIGURATION
############################

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cma16@iu.edu	# Add email here!
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=9:00:00
#SBATCH -A r00235

#SBATCH --job-name=aligned1_GSE196323	# CHANGE dataset ID
#SBATCH --output=/N/slate/cma16/logs/%x_%j.out	# CHANGE user cma16 here!
#SBATCH --error=/N/slate/cma16/logs/%x_%j.err	# CHANGE user cma16 here!

# Exit if any command fails, undefined variable is used, or pipe fails
set -euo pipefail

### -------------------- STEP 4 -------------------- ###
# load STAR to align trimmed RNA-seq reads to the reference genome
module load star

# Define directories and dataset ID
DATASET=GSE196323	# CHANGE dataset ID
SLATE_BASE=/N/slate/cma16	# CHANGE user cma16 here!
FASTQ_DIR=$SLATE_BASE/fastq/$DATASET # FASTQ directory
ALIGN_DIR=$SLATE_BASE/aligned1/$DATASET # alignment directory
TMP=/N/scratch/cma16 # Temporary folder

# Reference genome obtained from Ensembl (release 109 - assembly GRCh38)
GENOME_DIR=$SLATE_BASE/star_index #Human genome (reference)

# Create the alignment and temporary directory they do not exist
mkdir -p "$ALIGN_DIR"
mkdir -p "$TMP/$DATASET"

# List of samples
# Define single-end samples (one FASTQ file per sample)
SINGLE="SRR17933330 SRR17933331 SRR17933332 SRR17933333 SRR17933334 SRR17933335 SRR17933336 SRR17933337"

# Loop for STAR alignment (Single-end!)
# Align each single-end sample to the reference genome
for s in $SINGLE
do
	echo "Aligning ${s} ..."

	# Run STAR using a single FASTQ file (single-end reads)
	# Output includes a sorted BAM file and gene read counts
	STAR \
		--runThreadN $SLURM_CPUS_PER_TASK \
		--genomeDir "$GENOME_DIR" \
		--readFilesIn "$FASTQ_DIR/${s}_trimmed.fastq.gz" \
		--readFilesCommand zcat \
		--outFileNamePrefix "$ALIGN_DIR/${s}_" \
		--outTmpDir "$TMP/$DATASET/${s}_tmp_${SLURM_JOB_ID}" \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode GeneCounts

	# Confirm completion of alignment for the sample.
	echo "Finished ${s}"
done