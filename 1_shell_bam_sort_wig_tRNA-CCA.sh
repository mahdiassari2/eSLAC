#!/bin/sh
ref="hg38_tRNA_HC_intro_removed_CCA_upper_case"
refsuffix=".txt"
project="Amin_March2025_Nextseq"
reference=/project/taopan/mahdi/ref_seq/$ref$refsuffix
rm $reference.fai
echo $reference
index_file="${reference##*/}"
index_file="${index_file%%.[^.]*}"
echo $index_file

mkdir -p /project/taopan/mahdi/4_bam_sort_wig/$project/$index_file
mkdir -p /project/taopan/mahdi/5_tsv/$project/$index_file

for fullpath in /project/taopan/mahdi/3_bowtie2/$project/$index_file/*.sam
#for fullpath in /project/taopan/mahdi/4_bam_sort_wig/$project/*.sam
do
sleep 0.1
filename="${fullpath##*/}"    # Strip longest match of */ from start
#base="${filename%%.[^.]*}"    #Strip everything after the first period
base="${filename%%.[^.sam]}"
echo $base

echo "#!/bin/bash

#SBATCH --account=pi-taopan
#SBATCH --job-name=bam_sort_wig
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --partition=caslake
#SBATCH --mem-per-cpu=4000
#SBATCH --output=/project/taopan/mahdi/4_bam_sort_wig/$project/bam_sort_wig_$base.out
#SBATCH --error=/project/taopan/mahdi/4_bam_sort_wig/$project/bam_sort_wig_$base.err

#sleep 1

#module load midway2; 
module load python
module load java
module load python/cpython-3.8.5
module load samtools/1.13

filename="${fullpath##*/}"
base="${filename%%.[^.]*}"
echo $base
samtools view -bS -o /project/taopan/mahdi/4_bam_sort_wig/$project/$index_file/$base.bam /project/taopan/mahdi/3_bowtie2/$project/$index_file/$filename
sleep 0.1
samtools sort /project/taopan/mahdi/4_bam_sort_wig/$project/$index_file/$base.bam -o /project/taopan/mahdi/4_bam_sort_wig/$project/$index_file/$base.sort.bam
sleep 0.1
/project/taopan/mahdi/tools/IGV_2.8.0/igvtools count -z 5 -w 1 -e 250 --bases /project/taopan/mahdi/4_bam_sort_wig/$project/$index_file/$base.sort.bam /project/taopan/mahdi/4_bam_sort_wig/$project/$index_file/$base.wig $reference
sleep 0.1
python /project/taopan/mahdi/tools/wig_to_tsv_low_mem.py -i /project/taopan/mahdi/4_bam_sort_wig/$project/$index_file/$base.wig -r $reference -o /project/taopan/mahdi/5_tsv/$project/$index_file/$base.tsv

echo $base

">/project/taopan/mahdi/sbatch/jobfile.sbatch
sbatch /project/taopan/mahdi/sbatch/jobfile.sbatch

done
