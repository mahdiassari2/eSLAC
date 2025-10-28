#!/bin/bash
# ===========================================================
# SAM Counting Script
# ===========================================================
# Description:
#   Converts mapped SAM files into abundance TSV files using
#   a Python script. Generates one SBATCH file per sample.
#
# Instructions:
#   1. Update PROJECT, REFERENCE, BOWTIE_OUTPUT_DIR, and COUNTER_OUTPUT_DIR.
#   2. Ensure Python environment and sam_counter.py are accessible.
# ===========================================================

# -----------------------------
# User-defined parameters
# -----------------------------
PROJECT="PROJECT_NAME"
REFERENCE="/path/to/reference_file"  # Reference FASTA
BOWTIE_OUTPUT_DIR="/path/to/bowtie2_output"
COUNTER_OUTPUT_DIR="/path/to/sam_counter_output"
SBATCH_DIR="/path/to/sbatch_scripts"

SAM_COUNTER_SCRIPT="/path/to/sam_counter.py"

mkdir -p "$COUNTER_OUTPUT_DIR/$PROJECT"
mkdir -p "$SBATCH_DIR"

# Optional suffix for output files
OUTPUT_SUFFIX="_abundance.tsv"

# -----------------------------
# Iterate over SAM files
# -----------------------------
for FULLPATH in "$BOWTIE_OUTPUT_DIR"/$PROJECT/*/*.sam; do
    sleep 0.01

    FILENAME="${FULLPATH##*/}"           # Extract filename
    BASE="${FILENAME%%.sam}"             # Remove .sam extension
    INDEX_FILE=$(basename $(dirname "$FULLPATH"))

    echo "Processing SAM file: $FILENAME"

    SBATCH_FILE="$SBATCH_DIR/sam_counter_${BASE}.sbatch"

    cat <<EOF > "$SBATCH_FILE"
#!/bin/bash
#SBATCH --job-name=sam_count_${BASE}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --partition=YOUR_PARTITION
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mem=40000
#SBATCH --output=$COUNTER_OUTPUT_DIR/$PROJECT/$INDEX_FILE/sam_counter_${BASE}.out
#SBATCH --error=$COUNTER_OUTPUT_DIR/$PROJECT/$INDEX_FILE/sam_counter_${BASE}.err

# Load environment
module unload python
module load python

# Run SAM counting
python $SAM_COUNTER_SCRIPT -i $FULLPATH -o $COUNTER_OUTPUT_DIR/$PROJECT/$INDEX_FILE/$BASE$OUTPUT_SUFFIX

# Optional: remove mRNA (commented out)
# python3 /path/to/removemRNA.py -i $COUNTER_OUTPUT_DIR/$PROJECT/$INDEX_FILE/$BASE$OUTPUT_SUFFIX

EOF

    # Submit the SBATCH script
    sbatch "$SBATCH_FILE"
done
