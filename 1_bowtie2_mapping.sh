#!/bin/bash
# ===========================================================
# Bowtie2 Mapping and SAM Splitting Script
# Generic version for sharing
# ===========================================================
# Description:
#   Iterates over barcode directories containing FASTQ files,
#   aligns reads using Bowtie2, and splits SAM files for downstream analysis.
#
# Instructions:
#   1. Update paths at the top for your environment and reference files.
#   2. Ensure Bowtie2 and Python scripts are accessible.
#   3. Submit this script to generate SBATCH jobs (or run locally if desired).
# ===========================================================

# -----------------------------
# User-defined parameters
# -----------------------------

PROJECT="PROJECT_NAME"

REFERENCE="/path/to/bowtie2_index"

BARCODE_BASE_DIR="/path/to/barcode_dirs"

OUTPUT_BASE_DIR="/path/to/output_directory"

SBATCH_DIR="/path/to/sbatch_scripts"

SAM_SUFFIX=".sam"

BOWTIE2_BIN="/path/to/bowtie2"

SAM_SPLIT_SCRIPT="/path/to/sam_bin_split.py"

THREADS=8

# -----------------------------
# Prepare directories
# -----------------------------
mkdir -p "$OUTPUT_BASE_DIR"
mkdir -p "$SBATCH_DIR"

# -----------------------------
# Iterate over FASTQ files
# -----------------------------
for FULLPATH in "$BARCODE_BASE_DIR"/$PROJECT/*/*_2.fastq.gz; do
  sleep 0.1  # optional, to prevent filesystem overload

  # Extract sample info
  SAMPLE_DIR="${FULLPATH#*/barcode_dirs/}"
  SAMPLE_DIR="${SAMPLE_DIR%%/[^/]*}"
  BAR_DIR="${FULLPATH%/*.*}/"
  FILENAME="${FULLPATH##*/}"
  BASE="${FILENAME%%.[^.]*}"
  SAMPLE="${BASE%%_2*}"

  UNDERSCORE="_"
  echo "Processing sample: $SAMPLE_DIR$UNDERSCORE$SAMPLE"
  echo

  # Detect read orientation
  READ="${BAR_DIR#*/}"
  READ="${READ%%/[^/]*}"
  READ="${READ##*_}"

  if [ "$READ" == "read1" ]; then
      IN1="_1.fastq.gz"
      IN2="_2.fastq.gz"
  else
      IN1="_2.fastq.gz"
      IN2="_1.fastq.gz"
  fi

  # -----------------------------
  # Generate SBATCH script for this sample
  # -----------------------------
  SBATCH_FILE="$SBATCH_DIR/${SAMPLE_DIR}_${SAMPLE}.sbatch"
  
  cat <<EOF > "$SBATCH_FILE"
#!/bin/bash
#SBATCH --job-name=map_${SAMPLE}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=16000
#SBATCH --time=12:00:00
#SBATCH --output=$OUTPUT_BASE_DIR/$SAMPLE_DIR$UNDERSCORE$SAMPLE.out
#SBATCH --error=$OUTPUT_BASE_DIR/$SAMPLE_DIR$UNDERSCORE$SAMPLE.err

# Load environment
module load python  # optional
# Activate conda environment if needed
# source /path/to/conda.sh
# conda activate YOUR_ENV

# Bowtie2 alignment
$BOWTIE2_BIN -x $REFERENCE \\
  -U $BAR_DIR$SAMPLE$IN1 \\
  -S $OUTPUT_BASE_DIR/$SAMPLE_DIR$UNDERSCORE$SAMPLE$SAM_SUFFIX \\
  -q -p $THREADS --local --no-unal

# Split SAM into bins
python $SAM_SPLIT_SCRIPT \\
  -i $OUTPUT_BASE_DIR/$SAMPLE_DIR$UNDERSCORE$SAMPLE$SAM_SUFFIX \\
  -o $OUTPUT_BASE_DIR/ \\
  -breaks 0,10,20,30,40,50,60

echo "Finished sample: $SAMPLE_DIR$UNDERSCORE$SAMPLE"
EOF

  echo "Generated SBATCH script: $SBATCH_FILE"
  # Uncomment to submit
  # sbatch "$SBATCH_FILE"
done
