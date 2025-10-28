#!/bin/bash
# ===========================================================
# BAM Sorting, Wig Generation, and TSV Conversion Script
# Generic version for sharing
# ===========================================================
# Description:
#   Iterates over SAM files, converts to BAM, sorts, generates
#   coverage wig files with IGVTools, and converts wig to TSV.
#
# Instructions:
#   1. Update REFERENCE, PROJECT, and base directories.
#   2. Ensure samtools, java, IGVTools, and Python scripts are accessible.
#   3. This script will create SBATCH scripts for each sample.
# ===========================================================

# -----------------------------
# User-defined parameters
# -----------------------------

PROJECT="PROJECT_NAME"
REFERENCE="/path/to/reference_file"      # Reference fasta
BOWTIE_OUTPUT_DIR="/path/to/bowtie2_output"
BAM_OUTPUT_DIR="/path/to/bam_sort_wig"
TSV_OUTPUT_DIR="/path/to/tsv_output"
SBATCH_DIR="/path/to/sbatch_scripts"

IGVTOOLS_BIN="/path/to/igvtools"
WIG2TSV_SCRIPT="/path/to/wig_to_tsv_low_mem.py"

# Create output directories
mkdir -p "$BAM_OUTPUT_DIR/$PROJECT"
mkdir -p "$TSV_OUTPUT_DIR/$PROJECT"
mkdir -p "$SBATCH_DIR"

# -----------------------------
# Iterate over SAM files
# -----------------------------
for FULLPATH in "$BOWTIE_OUTPUT_DIR"/$PROJECT/*/*.sam; do
    sleep 0.1

    FILENAME="${FULLPATH##*/}"          # Extract filename
    BASE="${FILENAME%%.sam}"            # Remove .sam extension
    INDEX_FILE=$(basename $(dirname "$FULLPATH"))

    echo "Processing SAM file: $FILENAME"

    SBATCH_FILE="$SBATCH_DIR/bam_sort_wig_${BASE}.sbatch"

    cat <<EOF > "$SBATCH_FILE"
#!/bin/bash
#SBATCH --job-name=bam_sort_wig_${BASE}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --partition=YOUR_PARTITION
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mem-per-cpu=4000
#SBATCH --output=$BAM_OUTPUT_DIR/$PROJECT/bam_sort_wig_${BASE}.out
#SBATCH --error=$BAM_OUTPUT_DIR/$PROJECT/bam_sort_wig_${BASE}.err

# Load modules
module load python
module load java
module load samtools

# Convert SAM to BAM
samtools view -bS -o $BAM_OUTPUT_DIR/$PROJECT/$INDEX_FILE/$BASE.bam $FULLPATH
sleep 0.1

# Sort BAM
samtools sort $BAM_OUTPUT_DIR/$PROJECT/$INDEX_FILE/$BASE.bam -o $BAM_OUTPUT_DIR/$PROJECT/$INDEX_FILE/$BASE.sort.bam
sleep 0.1

# Generate coverage WIG
$IGVTOOLS_BIN count -z 5 -w 1 -e 250 --bases $BAM_OUTPUT_DIR/$PROJECT/$INDEX_FILE/$BASE.sort.bam $BAM_OUTPUT_DIR/$PROJECT/$INDEX_FILE/$BASE.wig $REFERENCE
sleep 0.1

# Convert WIG to TSV
python $WIG2TSV_SCRIPT -i $BAM_OUTPUT_DIR/$PROJECT/$INDEX_FILE/$BASE.wig -r $REFERENCE -o $TSV_OUTPUT_DIR/$PROJECT/$INDEX_FILE/$BASE.tsv

echo "Finished processing: $BASE"
EOF

    # Submit the SBATCH script
    sbatch "$SBATCH_FILE"
done
