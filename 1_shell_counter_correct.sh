#!/bin/sh
ref="mm10_tRNA202502"
refsuffic=".fasta"
project="miseq_20250721"

reference=/project/taopan/mahdi/ref_seq/$ref$refsuffix
rm $reference.fa

index_file="${reference##*/}"
index_file="${index_file%%.[^.]*}"
echo $index_file
mkdir -p /project/taopan/mahdi/6_sam_counter/$project/$index_file
suffix=_abundance.tsv

for fullpath in /project/taopan/mahdi/3_bowtie2/$project/$index_file/*.sam
do
sleep 0.01

filename="${fullpath##*/}"    # Strip longest match of */ from start
base="${filename%%.[^.]*}"    #Strip everything after the first period

echo "#!/bin/bash
#SBATCH --account=pi-taopan
#SBATCH --job-name=kallisto
#SBATCH --output=/project/taopan/mahdi/6_sam_counter/$project/$index_file/kallisto.out
#SBATCH --error=/project/taopan/mahdi/6_sam_counter/$project/$index_file/kallisto.err
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=40000

module unload python
#module load midway2
module load python


python /project/taopan/mahdi/tools/sam_counter.py -i $fullpath -o /project/taopan/mahdi/6_sam_counter/$project/$index_file/$base.tsv

sleep 1

#python3 ./removemRNA.py -i ./6_sam_counter/$index_file/$base.tsv

"> /project/taopan/mahdi/sbatch/jobfile.sbatch
sbatch /project/taopan/mahdi/sbatch/jobfile.sbatch

done
