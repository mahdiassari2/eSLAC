#!/bin/bash
ref="hg38_tRNA_HC_intro_removed_CCA_upper_case"
project="Amin_Nextseq_May2025"
reference=/project/taopan/mahdi/ref_seq/bowtie2_index/$ref/$ref

index_file="${reference##*/}"
mkdir -p /project/taopan/mahdi/3_bowtie2/$project/"$index_file"
suffix=.sam

for fullpath in /project/taopan/mahdi/0_barcode*/$project/*/*_2.txt.gz ; do
  sleep 0.1

  sample_dir="${fullpath#*/0_barcode*/$project/}" # gets the file and enclosing directory
  sample_dir="${sample_dir%%/[^/]*}" # remove the filename to get just the directory
  bar_directory="${fullpath%/*.*}/"

  filename="${fullpath##*/}"    # Strip longest match of */ from the start
  base="${filename%%.[^.]*}"    # Strip everything after the first period
  sample="${base%%_2*}"

  read="${bar_directory#*/}"
  read="${read%%/[^/]*}"
  read="${read##*_}"

  underscore="_"

#echo $bar_directory
#echo $sample_dir
#echo $base
#echo $sample
#echo $read
echo $sample_dir$underscore$sample
echo 

#Detect which read the barcode is on, and merge accordingly to orient the read as sense
if [ $read == "read1" ]
then
in1=_1.txt.gz
in2=_2.txt.gz
elif [ $read == "read2" ]
then
in1=_2.txt.gz
in2=_1.txt.gz
else
in1=_2.txt.gz
in2=_1.txt.gz
fi

echo "#!/bin/bash
#SBATCH --job-name=map
#SBATCH --account=pi-taopan
#SBATCH --partition=caslake
#SBATCH -o /project/taopan/mahdi/3_bowtie2/$project/$index_file/$sample_dir$underscore$sample.out
#SBATCH -e /project/taopan/mahdi/3_bowtie2/$project/$index_file/$sample_dir$underscore$sample.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=20000

sleep 1
#module load midway2
module unload python
module load python

/project/taopan/mahdi/tools/bowtie2-2.3.3.1-linux-x86_64/bowtie2 -x $reference -U $bar_directory$sample$in1 -S /project/taopan/mahdi/3_bowtie2/$project/$index_file/$sample_dir$underscore$sample$suffix -q -p 10 --local --no-unal
sleep 1
python /project/taopan/mahdi/tools/sam_bin_split.py -i /project/taopan/mahdi/3_bowtie2/$project/$index_file/$sample_dir$underscore$sample$suffix -o /project/taopan/mahdi/3_bowtie2/$project/$index_file/  -breaks 0,10,20,30,40,50,60

echo -e \"$sample_dir$underscore$sample\"
echo -e \"$bar_directory$sample$in1\"
echo -e \"\"
" >/project/taopan/mahdi/sbatch/jobfile.sbatch

  sbatch /project/taopan/mahdi/sbatch/jobfile.sbatch
done

