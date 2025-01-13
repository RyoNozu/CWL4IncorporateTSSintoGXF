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
    - BSgenome_data_package_seed_file 用にゲノムfastaファイルを sequence ID ごとに分割  
    - seqkit を使用  
        - Version: 2.8.2  
    - option setting  
        - -i ; sequence id で分割  
        - -o ; output directory の指定  

- 