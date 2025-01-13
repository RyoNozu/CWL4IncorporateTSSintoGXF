#! /bin/zsh

# This script is used to map the reads to the Htri genome using STAR
# check ulimit "ulimit -a" 
# before running # set "ulimit -n 2048"

sample=($(ls *_trim_R1.fq.gz | sed 's/_trim_R1.fq.gz//'))

ulimit -n 2048

echo "run STAR for Htri selected reference genome and gtf file" `date '+%Y/%m/%d %H:%M:%S'` >> run_star.log
STAR --version >> run_star.log 2>&1

for i in ${sample[@]};do
STAR --genomeDir ./${index_dir}/ \
--readFilesCommand zless \
--readFilesIn ${i}_trim_R1.fq.gz ${i}_trim_R2.fq.gz \
--runThreadN `sysctl -n hw.physicalcpu_max` \
--alignEndsType Extend5pOfRead1 \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${i} >> run_star.log 2>&1
done
