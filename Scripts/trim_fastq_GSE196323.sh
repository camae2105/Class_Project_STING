#!/bin/bash

############################
# SLURM JOB CONFIGURATION
############################

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cma16@iu.edu	# Add your email here!
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH -A r00235

# Set the job name and output/error log files for this trimming step
#SBATCH --job-name=trimFastq_GSE196323	# CHANGE dataset ID
#SBATCH --output=/N/slate/cma16/logs/%x_%j.out  # CHANGE user
#SBATCH --error=/N/slate/cma16/logs/%x_%j.err	# CHANGE user

### -------------------- STEP 3 -------------------- ###
# Trim fastq files in SLATE folder

# Load the fastp module for read trimming and quality control
module load fastp

# Define folders
# Set the dataset ID and define the FASTQ input directory
DATASET=GSE196323	# CHANGE depending on your dataset ID
SLATE_BASE=/N/slate/cma16	# CHANGE the user cma16 here
FASTQ_DIR=$SLATE_BASE/fastq/$DATASET

# List of samples
# Separate paired-end and single-end samples before trimming
PAIRED="SRR17933314 SRR17933315 SRR17933316 SRR17933317 SRR17933318 SRR17933319 SRR17933320 SRR17933321 SRR17933322 SRR17933323 SRR17933324 SRR17933325 SRR17933326 SRR17933327 SRR17933328 SRR17933329"
SINGLE="SRR17933330 SRR17933331 SRR17933332 SRR17933333 SRR17933334 SRR17933335 SRR17933336 SRR17933337"

###### For paired-end reads ######

# Trim paired-end FASTQ files and generate quality-filtered outputs
echo "===== TRIMMING READS (PAIRED-END, EXISTING FASTQ.GZ) ====="

for s in $PAIRED
do
	echo "Trimming $s"

	# Run fastp on paired-end reads and generate QC reports.
	fastp \
		-i "$FASTQ_DIR/${s}_1.fastq.gz" \
		-I "$FASTQ_DIR/${s}_2.fastq.gz" \
		-o "$FASTQ_DIR/${s}_1_trimmed.fastq.gz" \
		-O "$FASTQ_DIR/${s}_2_trimmed.fastq.gz" \
		-w $SLURM_CPUS_PER_TASK \
		--detect_adapter_for_pe \
		--html "$FASTQ_DIR/${s}_fastp.html" \
		--json "$FASTQ_DIR/${s}_fastp.json"

	# Check that both trimmed files were created before deleting the original files
	if [[ -s "$FASTQ_DIR/${s}_1_trimmed.fastq.gz" && -s "$FASTQ_DIR/${s}_2_trimmed.fastq.gz" ]]; then
		rm "$FASTQ_DIR/${s}_1.fastq.gz"
		rm "$FASTQ_DIR/${s}_2.fastq.gz"
	else
		echo "ERROR: Trimming failed for $s"
		exit 1
	fi
done

###### For single-end reads ######

# Trim single-end FASTQ files and generate quality-filtered outputs
echo "===== TRIMMING READS (SINGLE-END, EXISTING FASTQ.GZ) ====="

for s in $SINGLE
do
	echo "Trimming $s"

	# Run fastp on single-end reads and generate QC reports.
	fastp \
		-i "$FASTQ_DIR/${s}.fastq.gz" \
		-o "$FASTQ_DIR/${s}_trimmed.fastq.gz" \
		-w $SLURM_CPUS_PER_TASK \
		--html "$FASTQ_DIR/${s}_fastp.html" \
		--json "$FASTQ_DIR/${s}_fastp.json"

	# Check that the trimmed file was created before deleting the original file
	if [[ -s "$FASTQ_DIR/${s}_trimmed.fastq.gz" ]]; then
		rm "$FASTQ_DIR/${s}.fastq.gz"
	else
		echo "ERROR: Trimming failed for $s"
		exit 1
	fi
done
