#!/bin/bash

# ------------------------------
# SLURM JOB CONFIGURATION
# ------------------------------
# Define job resources and logging
#SBATCH --job-name=counts_TPM_GSE196323
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cma16@iu.edu	# add your email here
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH -A r00235
#SBATCH --output=/N/slate/cma16/logs/%x_%j.out # CHANGE cma16 here
#SBATCH --error=/N/slate/cma16/logs/%x_%j.err	# CHANGE cma16 here

# Exit if any command fails, undefined variable is used, or pipe fails
set -euo pipefail

# Load required modules
# subread → for featureCounts (gene counting)
# python → to run downstream normalization scripts
module load subread
module load python

# ------------------------------
# DEFINE DIRECTORIES AND FILES
# ------------------------------

# Define dataset ID and directories
DATASET=GSE196323	# CHANGE dataset ID!
SLATE_BASE=/N/slate/cma16 # Base directory - CHANGE user cma16 here!!!
ALIGN_DIR=$SLATE_BASE/aligned1/$DATASET # With BAM files
COUNT_DIR=$SLATE_BASE/counts1/$DATASET # To store counts and expression matrices
REF_DIR=$SLATE_BASE/ref # Reference directory 
GTF=$REF_DIR/Homo_sapiens.GRCh38.109.gtf # GTF files for counting reads
GENE_MAP=$REF_DIR/gene_id_to_name.tsv # Mapping to add gene names

# Example of format for gene_id_to_name.tsv:
# ENSG00000000003	TSPAN6
# ENSG00000000005	TNMD
# ENSG00000000419	DPM1
# ENSG00000000457	SCYL3

# Create output directory if it does not exist
mkdir -p "$COUNT_DIR"

# ------------------------------
# DEFINE OUTPUT FILES
# ------------------------------

# Raw count files
PAIRED_COUNTS=$COUNT_DIR/paired_counts.txt
SINGLE_COUNTS=$COUNT_DIR/single_counts.txt

# Combined counts file
MERGED_COUNTS=$COUNT_DIR/all_counts_${DATASET}.txt
CLEAN_COUNTS=$COUNT_DIR/all_counts_${DATASET}_cleaned.txt

# Expression matrices
TPM_OUT=$COUNT_DIR/all_TPM_${DATASET}.tsv
LOG2TPM_OUT=$COUNT_DIR/all_log2TPM_${DATASET}.tsv

# Expression matrices with gene names
TPM_NAMED=$COUNT_DIR/all_TPM_${DATASET}_withNames.tsv
LOG2TPM_NAMED=$COUNT_DIR/all_log2TPM_${DATASET}_withNames.tsv

# ------------------------------
# SAMPLE LISTS
# ------------------------------

# Paired-end samples
PAIRED="SRR17933298 SRR17933299 SRR17933300 SRR17933301 SRR17933302 SRR17933303 SRR17933304 SRR17933305 SRR17933306 SRR17933307 SRR17933308 SRR17933309 SRR17933310 SRR17933311 SRR17933312 SRR17933313 SRR17933314 SRR17933315 SRR17933316 SRR17933317 SRR17933318 SRR17933319 SRR17933320 SRR17933321 SRR17933322 SRR17933323 SRR17933324 SRR17933325 SRR17933326 SRR17933327 SRR17933328 SRR17933329"

# Single-end samples
SINGLE="SRR17933330 SRR17933331 SRR17933332 SRR17933333 SRR17933334 SRR17933335 SRR17933336 SRR17933337"

# ------------------------------
# STEP 1: GENE COUNTING
# ------------------------------

# Run featureCounts for paired-end samples
echo "Running featureCounts for paired-end samples..."

# Add -p for paired-end mode
featureCounts \
	-T $SLURM_CPUS_PER_TASK \ 
	-p \
	-t exon \
	-g gene_id \
	-a "$GTF" \
	-o "$PAIRED_COUNTS" \
	$(for s in $PAIRED; do echo "$ALIGN_DIR/${s}_Aligned.sortedByCoord.out.bam"; done)

# Run featureCounts for single-end samples
echo "Running featureCounts for single-end samples..."

# Remove -p for single-end mode
featureCounts \
	-T $SLURM_CPUS_PER_TASK \
	-t exon \
	-g gene_id \
	-a "$GTF" \
	-o "$SINGLE_COUNTS" \
	$(for s in $SINGLE; do echo "$ALIGN_DIR/${s}_Aligned.sortedByCoord.out.bam"; done)

# ------------------------------
# STEP 2: MERGE COUNT FILES
# ------------------------------

echo "Merging paired and single count tables..."

# Extract metadata columns (Geneid, Chr, etc.)
cut -f1-6 "$PAIRED_COUNTS" > "$MERGED_COUNTS"

# Merge sample columns from paired and single counts
paste "$MERGED_COUNTS" \
	<(cut -f7- "$PAIRED_COUNTS") \
	<(cut -f7- "$SINGLE_COUNTS") \
	> "$COUNT_DIR/tmp.txt"

# Replace merged file
mv "$COUNT_DIR/tmp.txt" "$MERGED_COUNTS"

# ------------------------------
# STEP 3: CLEAN COLUMN NAMES
# ------------------------------

echo "Cleaning headers..."

# Remove file paths and simplify sample names
awk 'NR==1 { print; next }
NR==2 {
	for (i=1; i<=NF; i++) {
		gsub(".*/","",$i)
		gsub("_Aligned.sortedByCoord.out.bam","",$i)
	}
	print
	next
}
NR>2 { print }' OFS="\t" "$MERGED_COUNTS" > "$CLEAN_COUNTS"

# ------------------------------
# STEP 4: NORMALIZATION (CPM)
# ------------------------------

echo "Computing CPM..."

# Convert raw counts to Counts Per Million (CPM)
python /N/slate/cma16/python_scripts/getCPMmatrix.py \
	"$CLEAN_COUNTS" \
	"$COUNT_DIR/all_CPM_${DATASET}.tsv"

# ------------------------------
# STEP 5: NORMALIZATION (TPM)
# ------------------------------

echo "Computing TPM..."

# Compute TPM and log2(TPM+1) matrices
python /N/slate/cma16/python_scripts/twoTPM_matrix.py \
	"$CLEAN_COUNTS" \
	"$TPM_OUT" \
	"$LOG2TPM_OUT"

# ------------------------------
# STEP 6: ADD GENE NAMES
# ------------------------------

echo "Adding gene names..."

# Map gene IDs to gene symbols for easier interpretation
python /N/slate/cma16/python_scripts/add_gene_names.py \
	"$TPM_OUT" \
	"$GENE_MAP" \
	"$TPM_NAMED"

python /N/slate/cma16/python_scripts/add_gene_names.py \
	"$LOG2TPM_OUT" \
	"$GENE_MAP" \
	"$LOG2TPM_NAMED"

# ------------------------------
# FINAL MESSAGE
# ------------------------------

echo "All steps completed successfully."