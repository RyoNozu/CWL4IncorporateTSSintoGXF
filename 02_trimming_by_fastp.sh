#! /bin/zsh

# run fastp for triming only adapter

sample=($(ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//'))
# -q: quality threshold
# -Q: If this option is specified, quality filtering is disabled
# -l: length threshold
# -z: compression level for gzip output. 9 is smallest
# -w: number of threads

for i in ${sample[@]};
do
fastp -i ${i}_R1.fastq.gz -I ${i}_R2.fastq.gz \
-o ${i}_trim_R1.fq.gz -O ${i}_trim_R2.fq.gz \
-h ${i}_trim_report.html -j ${i}_trim_report.json \
-Q -z 9 -w 12 >>run_fastp.log 2>&1
done
