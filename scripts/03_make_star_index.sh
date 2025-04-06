#! /bin/zsh

# this script is used to mapping the reads to Htri reference genome
#
# 引数の指定例
# $1=star_genome_idx #出力用のディレクトリ名
# $2="Halichoeres_trimaculatus-hifiasm-3ddna-v1.1.edit.sfmasked.selected.fa" # reference ゲノム配列（fasta）
# $3="braker_correctID_v3.gtf" # reference ゲノムのアノテーション情報（gtf）
# 


echo "Start to make STAR index for Htri reference genome gtf file" `date '+%Y/%m/%d %H:%M:%S'` >> prepare_star-index.log
STAR --version >> prepare_star-index.log 2>&1
mkdir $1 #出力用のディレクトリを作成
gzip -d -k ./$2.gz
STAR --runMode genomeGenerate \
--genomeDir $1/ \
--genomeFastaFiles $2 \
--sjdbGTFfile $3 \
--sjdbOverhang 100 \
--runThreadN `sysctl -n hw.physicalcpu_max` >> prepare_star-index.log 2>&1
