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
# List the paired-end samples to be aligned in this job
PAIRED="SRR17933314 SRR17933315 SRR17933316 SRR17933317 SRR17933318 SRR17933319 SRR17933320 SRR17933321 SRR17933322 SRR17933323 SRR17933324 SRR17933325 SRR17933326 SRR17933327 SRR17933328 SRR17933329"

# Loop for STAR alignment (Paired-end!)
# Align each paired-end sample to the genome
for s in $PAIRED
do
	echo "Aligning ${s} ..."

	# Run STAR using the two trimmed FASTQ files for each paired-end sample
	# Output includes a sorted BAM file and gene read counts
	STAR \
		--runThreadN $SLURM_CPUS_PER_TASK \
		--genomeDir "$GENOME_DIR" \
		--readFilesIn "$FASTQ_DIR/${s}_1_trimmed.fastq.gz" "$FASTQ_DIR/${s}_2_trimmed.fastq.gz" \
		--readFilesCommand zcat \
		--outFileNamePrefix "$ALIGN_DIR/${s}_" \
		--outTmpDir "$TMP/$DATASET/${s}_tmp_${SLURM_JOB_ID}" \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode GeneCounts

	# Confirm that the sample finished aligning before moving to the next one
	echo "Finished ${s}"
done