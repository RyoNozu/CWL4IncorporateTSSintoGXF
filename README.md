# CWL4gtf_update

- private repository for creating CWL for updating gtf  

## 内容  

- Input files を取り込んで、TSS を決定し、TSS に基づいて annotation 5'UTR を gtf file に付加しアップデートした gtf file を output file として出力する  

### Input files  
- CAGE-seq Read  
- reference genome (fasta)  
- gtf file  
- BSgenome_data_package_seed_file (.txt)  

### Output file  
- updated gtf file (.gtf)  

***

#### scripts for each step  

- 01_split_genome_seqs_4_prepare_BSgenome_data_package_seed_file.sh  
  
      % zsh 01_split_genome_seqs_4_prepare_BSgenome_data_package_seed_file.sh "genome file"

    - BSgenome_data_package_seed_file 用にゲノムfastaファイルを sequence ID ごとに分割  
    - seqkit を使用  
        - Version: 2.8.2  
    - option setting  
        - -i ; sequence id で分割  
        - -o ; output directory の指定  


- 02_trimming_by_fastp.sh  
  
      % zsh 02_trimming_by_fastp.sh  

    - CAGE-Seq read を fastp でトリミング  
    - default  
    - adaptor のトリミングが主  
    - fastp  
        - Version: 0.23.4  
    - Input: cage-seq read (.fastq.gz)  

- 03_make_star_index.sh  
  
      % zsh 03_make_star_index.sh "出力用のディレクトリ名" "reference genome" "gtf file"  

    - マッピングソフト STAR の index を作成  
    - 引数1の値でディレクトリを作成し、そこに出力  
    - STAR  
        - Version:  
    - 

- 04_run_STAR_4_CAGE-Seq_analysis.sh  
  
      zsh 04_run_STAR_4_CAGE-Seq_analysis.sh "03_ の引数1"  

    - STAR によるマッピング  
    - sorted bam を出力  
    - 出力形式 (*Aligned.sortedByCoord.out.bam)  
    - option setting  
        - 