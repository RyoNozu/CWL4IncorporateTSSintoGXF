# CWL4updating_gtf

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
    - adaptor のトリミングが主目的  
    - fastp  
        - Version: 0.23.4  
    - Input: cage-seq read (.fastq.gz)  

- 03_make_star_index.sh  
  
      % zsh 03_make_star_index.sh "出力用のディレクトリ名" "reference genome" "gtf file"  

    - マッピングソフト STAR の index を作成  
    - 引数1の値でディレクトリを作成し、そこに出力  
    - STAR  
        - Version: 2.7.11b  
    - 

- 04_run_STAR_4_CAGE-Seq_analysis.sh  
  
      zsh 04_run_STAR_4_CAGE-Seq_analysis.sh "03_ の引数1"  

    - STAR によるマッピング  
    - sorted bam を出力  
    - 出力形式 (*Aligned.sortedByCoord.out.bam)  
    - option setting  
        - --alignEndsType Extend5pOfRead1  
          
              Note from TSSr README: If you ran STAR with the default alignment parameters (without using --alignEndsType Extend5pOfRead1), please ensure that you set softclippingAllowed = TRUE when running TSSr.  
              => Mapping時 (by STAR) ↑オプションを付けておく  
        -  

- 05_Install_and_load_packages_4_TSSr.R  

- 06_Generate_BSgenome_Package_4_target_species.r  

    - seed file の準備 (下記にHtriの1例)  
      
          Package: BSgenome.Htrimaculatus.inhouse.Htriv1  
          Title: Full genomic sequences for Halichoeres trimaculatus (Htriv1)  
          Description: Full genomic sequences for Halichoeres trimaculatus (inhouse version  
          Htriv1).  
          Version: 1.0.0  
          organism: Halichoeres trimaculatus  
          common_name: three-spot wrasse  
          genome: Htriv1  
          provider: inhouse  
          release_date: 2023/07/01  
          source_url: https://github.com/RyoNozu/Htrimaculatus_genome  
          organism_biocview: Halichoeres_trimaculatus  
          BSgenomeObjname: Htrimaculatus  
          circ_seqs: 0  
          seqnames: paste("HiC_scaffold_", c(1:227), sep="")  
          seqs_srcdir: /Volumes/project＿HtriTSS＿BU/Project_HtriTSS/ref_genome  
          info/seqs_by_scaffold  
          seqfiles_prefix: Halichoeres_trimaculatus-hifiasm-3ddna-v1.1.edit.part_  
          seqfiles_suffix: .fa.gz  
    
    - 染色体 (scaffold) ごとに分割した genome sequence が必要  
        - 01_split_xxx で対応済  

- 07_Import_required_files_and_TSS_calling.R  
