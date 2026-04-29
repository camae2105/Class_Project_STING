#!/bin/bash

############################
# SLURM JOB CONFIGURATION
############################

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cma16@iu.edu	# Add your email here
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=10:00:00
#SBATCH -A r00235

# Change dataset ID in next line for log files!!!
#SBATCH --job-name=getFastq_GSE196323	# CHANGE dataset ID
#SBATCH --output=/N/slate/cma16/logs/%x_%j.out	# CHANGE user cma16 here!
#SBATCH --error=/N/slate/cma16/logs/%x_%j.err	# CHANGE user cma16 here!

### -------------------- STEP 1 -------------------- ###
# Define directories and dataset

# Set the dataset ID for this run.
DATASET=GSE196323 #CHANGE this depending on your dataset ID!!!

# Define the scratch base directory used for temporary files
SCRATCH_BASE=/N/scratch/cma16/project #CHANGE the user cma16 here!

# Define the slate base directory used for long-term storage
SLATE_BASE=/N/slate/cma16	#CHANGE the user here!

# Define the temporary FASTQ directory in scratch space
FASTQ_DIR=$SCRATCH_BASE/fastq/$DATASET #Convert them in scratch

# Define the final FASTQ storage directory in slate.
NEW_FASTQ_DIR=$SLATE_BASE/fastq/$DATASET #Store them in slate

# Download .sra files via prefetch
# Load the SRA Toolkit module so prefetch and fasterq-dump are available
module load sra-toolkit/3.0.5 #Load this to use prefetch and fasterq-dump

# Check that fasterq-dump is available in the current environment
which fasterq-dump

# Print dataset information to make sure they're the right ones
echo "Dataset: $DATASET"
echo "FASTQ directory: $FASTQ_DIR"

Create the temporary and final FASTQ directory if they do not exist
mkdir -p "$FASTQ_DIR"
mkdir -p "$NEW_FASTQ_DIR"

# List all SRR accessions included in this dataset
SRR="SRR17933298 SRR17933299 SRR17933300 SRR17933301 SRR17933302 SRR17933303 SRR17933304 SRR17933305 SRR17933306 SRR17933307 SRR17933308 SRR17933309 SRR17933310 SRR17933311 SRR17933312 SRR17933313 SRR17933314 SRR17933315 SRR17933316 SRR17933317 SRR17933318 SRR17933319 SRR17933320 SRR17933321 SRR17933322 SRR17933323 SRR17933324 SRR17933325 SRR17933326 SRR17933327 SRR17933328 SRR17933329 SRR17933330 SRR17933331 SRR17933332 SRR17933333 SRR17933334 SRR17933335 SRR17933336 SRR17933337"

# Download each SRA run into the scratch directory
# This block is currently commented out because the SRA files are already available
# for s in $SRR
# do
#	echo "Downloading ${s}..."
#	prefetch --output-directory $FASTQ_DIR $s
# done
#
# Move the .sra files out of the SRR folders so they're all in the same folder
# mv $FASTQ_DIR/SRR*/SRR*.sra $FASTQ_DIR
#
# Remove the temporary subdirectories created by prefetch
# rmdir $FASTQ_DIR/SRR*

### -------------------- STEP 2 -------------------- ###
# Get FAST files from SRA files

# Separate runs into paired-end and single-end samples when necessary
# Check if they're paired-end or single-end before this step!
SINGLE="SRR17933330 SRR17933331 SRR17933332 SRR17933333 SRR17933334 SRR17933335 SRR17933336 SRR17933337"
PAIRED="SRR17933314 SRR17933315 SRR17933316 SRR17933317 SRR17933318 SRR17933319 SRR17933320 SRR17933321 SRR17933322 SRR17933323 SRR17933324 SRR17933325 SRR17933326 SRR17933327 SRR17933328 SRR17933329"

# Convert SRA → FASTQ for single-end reads
for s in $SINGLE
do
	echo "Converting ${s} (single-end)..."

	# Stop if the expected .sra file is missing
	if [[ ! -f "$FASTQ_DIR/${s}.sra" ]]; then
		echo "ERROR: $FASTQ_DIR/${s}.sra not found"
		exit 1
	fi

	# Convert the SRA file into a FASTQ file
	fasterq-dump \
		-e 8 \
		--temp "$FASTQ_DIR" \
		-O "$FASTQ_DIR" \
		$s

	# Compress the resulting FASTQ file
	pigz -p 8 "$FASTQ_DIR/${s}.fastq"

	# Move the compressed FASTQ file to long-term storage
	mv "$FASTQ_DIR/${s}.fastq.gz" "$NEW_FASTQ_DIR/"
#	rm "$FASTQ_DIR/${s}.sra" # Better remove them manually just in case
done

# Convert SRA → FASTQ for paired-end reads
for s in $PAIRED
do
	echo "Converting ${s} (paired-end)..."

	# Stop if the expected .sra file is missing
	if [[ ! -f "$FASTQ_DIR/${s}.sra" ]]; then
		echo "ERROR: $FASTQ_DIR/${s}.sra not found"
		exit 1
	fi

	# Convert the SRA file into paired FASTQ files
	fasterq-dump \
		--split-files \
		-e 8 \
		--temp "$FASTQ_DIR" \
		-O "$FASTQ_DIR" \
		$s

	# Compress the reads
	pigz -p 8 "$FASTQ_DIR/${s}_1.fastq"
	pigz -p 8 "$FASTQ_DIR/${s}_2.fastq"

	# Move the compressed FASTQ files to long-term storage
	mv "$FASTQ_DIR/${s}_1.fastq.gz" "$NEW_FASTQ_DIR/"
	mv "$FASTQ_DIR/${s}_2.fastq.gz" "$NEW_FASTQ_DIR/"
#	rm "$FASTQ_DIR/${s}.sra"
done