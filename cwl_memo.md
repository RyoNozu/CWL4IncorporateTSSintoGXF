# cwl-memo

## CWLを記述する手順

- CWLを書く時に以下の流れで記載をしています

### 1. zatsu-cwl-generatorでざっくりCWLファイルを出力する

```bash
zatsu-cwl-generator "seqkit ./Data/Halichoeres_trimaculatus/Halichoeres_trimaculatus-hifiasm-3ddna-v1.1.edit.fna.gz -i -O ./seqs_by_scaffold" > ./Tools/01_split_genome_seqs.cwl
```

### 2. 手動で修正，チェックする

- CWLの記法のバージョンを1.2にする
- 以下のコマンドでvalidateして書き方に変なところがないか確認する

```bash
cwltool --validate ./Tools/01_split_genome_seqs.cwl 
```

### 3. チェックが通ったら，ワークフローエンジンのcwltoolで実行する

```bash
cwltool --debug --cachedir ./cwl_cache/ --outdir ./Data/ ./Tools/01_split_genome_seqs.cwl
```

&nbsp;

&nbsp;

&nbsp;

## 2025/02/08 memo by yonezawa

- fastpのプロセスで複数ファイルの処理を行うCWLファイルを新たに作成
- サブワークフローを作成して，複数ファイルの処理を行う
- YAMLファイルにファイルのパスを記載して，サブワークフローを実行


```bash
cwltool --debug --outdir ./out/ ./workflow/01_trimming_fastq_subworkflow.cwl ./config/01_trimming_fastq_files.yml
```

&nbsp;

## 2025/02/09 memo by yonezawa

- STARでindexを作成するファイルを作成
- reference genomeのfastaファイルをdecompressする処理とSTARでindexを作成する処理を分離

```bash
cwltool --debug --outdir ./Data/Halichoeres_trimaculatus/ ./Tools/03_pigz.cwl
cwltool --debug --outdir ./out/ --cachedir ./cwl_cache/ ./Tools/04_make_star_index.cwl
```

&nbsp;

- STARによるマッピングプロセスも `ScatterFeatureRequirement` を使ってサブワークフローにすることで複数ファイルに対応
- なぜか `--readFilesCommand` のオプションを `zless` にするとエラーが出るので，`--readFilesCommand zcat` にしている

```bash
cwltool --debug --outdir ./out/ --cachedir ./cwl_cache/ ./workflow/02_star4cageseq_analysis_subworkflow.cwl ./config/02_star4cageseq_analysis.yml
```

&nbsp;

## 2025/03/13 memo by NZR

オプション指定での実行でエラー

```bash
cwltool --debug --outdir ./out/ --cachedir ./cwl_cache/ ~/Google\ Drive/その他の 
  ソコン/マイ\ MacBook\ Pro/Google_drive/Project_non-codingRegionAnnotation/CWL4gtf_update/CWL4gtf_update/Tools/04_make_star_index.cwl --reference_genome ./Data/Danio_rerio/GCF_000002035.6_GRCz11_genomic.fna --reference_genome_annotation ./Data/Danio_rerio/GCF_000002035.6_GRCz11_genomic.gtf
```

```DEBUG & ERROR message
DEBUG Parsed job order from command line: {
    "__id": "/Users/rnz/Google Drive/\u305d\u306e\u4ed6\u306e\u30d1\u30bd\u30b3\u30f3/\u30de\u30a4 MacBook Pro/Google_drive/Project_non-codingRegionAnnotation/CWL4gtf_update/CWL4gtf_update/Tools/04_make_star_index.cwl",
    "output_dir_name": "star_genome_idx",
    "reference_genome": {
        "class": "File",
        "location": "file:///Volumes/externalSSD4M3Max/test_dir_4_update-gtf_cwl/Data/Danio_rerio/GCF_000002035.6_GRCz11_genomic.fna"
    },
    "reference_genome_annotation": {
        "class": "File",
        "location": "file:///Volumes/externalSSD4M3Max/test_dir_4_update-gtf_cwl/Data/Danio_rerio/GCF_000002035.6_GRCz11_genomic.gtf"
    },
    "star_threads": 16
}
ERROR Workflow error:
Expected value of 'reference_genome' to have format 'http://edamontology.org/format_1929' but
 File has no 'format' defined: {
    "class": "File",
    "location": "file:///Volumes/externalSSD4M3Max/test_dir_4_update-gtf_cwl/Data/Danio_rerio/GCF_000002035.6_GRCz11_genomic.fna",
    "size": 1700419557,
    "basename": "GCF_000002035.6_GRCz11_genomic.fna",
    "nameroot": "GCF_000002035.6_GRCz11_genomic",
```

&nbsp;

## 2025/03/19 memo by yonezawa

### Workflow 構築基本方針

- single-end と pair-end に対応した処理を実行することに決定
- パーツである `CommandLineTool`定義のファイルをそれぞれ SEとPEで作成
- サブのパーツである `subworkflow` のファイルもそれぞれSE、PEで作成
- (悩み中) すべてのプロセスを含む `workflow`定義のファイルに `when` 式を組み込むか、それともこちらもSEとPEのバージョンを両方作るか...

&nbsp;

## 2025/03/20 memo by yonezawa

- `./Tools` ディレクトリにある `CommandLineTool`定義のファイル (ワークフローのパーツとなるもの) のうち、`fastp`と`STAR`の処理はそれぞれsingle-endに対応させた
- `cwltool --validate` を実行して一応OKだったものの、まだ実際のファイルで実行していない
- `./workflow`ディレクトリにある `subworkflow`という名前がついたファイルについてもsingle-endに対応させたものを使用
- こちらについてもまだ実際のファイルで実行していない

&nbsp;

### `workflow`定義のファイルについて

- 当初、`when`式を組み込んだワークフローのファイルを作成しようとかんがえていたが、以下のような懸念点がある
    1. CWLではシンプルなワークフローを構築することを推奨している
    2. この機能はすべてのワークフローエンジンには搭載されていないため、再現性、互換性に難あり
    3. シンプルな条件分岐なのに対して、複雑な対応が迫られる

- 以上の点からペアエンドとシングルエンドに対応したファイルをそれぞれ作成する
- READMEとかには ペアエンド用はこれ、シングルエンド用はこれ、というふうな指示を書いておく 


- seqkit -> STAR mapping までの処理 (paired-end) は一応成功

```bash
cwltool --debug --outdir ./test/ --cachedir ./cwl_cache/ ./workflow/cageseq_gtf_update_pe.cwl ./config/Workflow_config/cageseq_gtf_update_pe.yml
```

&nbsp;

## 2025/04/12 memo by yonezawa

### 全体の流れについて

- TSSrから09-02までの処理はこちらの環境で実行を確認し、野津さんにも確認いただいた
- 各スクリプトはもしかしたらこちらで少し改良してCWL化するかも
- 


```bash
Rscript 05--08-combined_exec_TSSr.R ./Data/Halichoeres_trimaculatus/braker_correctID_v3.gtf 2 1 1 A B

python3 09-01_join_all_assignedClusters.py *.assignedClusters.txt # ファイルを一つずつ指定することもできる

python3 09-02_extract_tss-feature_from_all-joined.assignedClusters_then_uniq.tss-feature.py 

# これはまだ改良するかもとのこと (by 野津さん)
# 最終的に野津さんが修正し、OKになった
python3 10_draft_update_gtf_add_all.assignedClusters_TSS-feature_v20250403-01.py Data/Halichoeres_trimaculatus/braker_correctID_v3.gtf all_tss_feature_uniq.gene.tsv test.gtf > log.txt
```

### GTFとGFFの分岐について

- これもできたらやる
- ワークフローとして別々に作ってもいいかも

&nbsp;

&nbsp;

### Rscriptの修正

- CWLでも実行できるように`optparse` を使用した方法に変更
- BAMファイルは `,`ですべてのファイルを渡すようにすることで対応
- おそらくCWLでも `itemSeparator` を使用してarrayになっているファイル一覧をつなげられると考えられる
- [./scripts/05-08_combined_exec_TSSr_modified.R](./scripts/05-08_combined_exec_TSSr_modified.R) に変更

```bash
Rscript ./scripts/05-08_combined_exec_TSSr_modified.R \
--refSource ./Data/Halichoeres_trimaculatus/braker_correctID_v3.gtf \
--seedFile ./Data/Halichoeres_trimaculatus/BSgenome.Htrimaculatus.Htriv1.1-seed_2.txt \
--group_count 2 \
--group_sizes 1,1 \
--prefixes A,B \
--threads 16 \
--inputFilesType bamPairedEnd \
--organismName Halichoeres trimaculatus \
--bam_files ./MK.F1_R1_trim.fq_Aligned.sortedByCoord.out.bam,./MK.F1.Mix_R1_trim.fq_Aligned.sortedByCoord.out.bam 
```

&nbsp;

&nbsp;

## 20250413

- RscriptをCWLで実行することを見据えて改造
- metadataを与えてグループにするBEDファイルを明示的に指定

```bash
Rscript ./scripts/05-08_combined_exec_TSSr_modified.R \
--refSource ./Data/Halichoeres_trimaculatus/braker_correctID_v3.gtf \
--seedFile ./Data/Halichoeres_trimaculatus/BSgenome.Htrimaculatus.Htriv1.1-seed_2.txt \
--threads 16 \
--inputFilesType bamPairedEnd \
--organismName Halichoeres trimaculatus \
--bam_dir ./out \
--metadata ./Data/Halichoeres_trimaculatus/sample_metadata.csv
```

&nbsp;

- 途中まで実行してみてBAMファイルを直接指定するバージョンでも実行できそうなので、途中で中断し、このパターンでCWL化に進む

```bash
Rscript ./scripts/05-08_combined_exec_TSSr_modified.R \
--refSource ./Data/Halichoeres_trimaculatus/braker_correctID_v3.gtf \
--seedFile ./Data/Halichoeres_trimaculatus/BSgenome.Htrimaculatus.Htriv1.1-seed_2.txt \
--threads 16 \
--inputFilesType bamPairedEnd \
--organismName Halichoeres trimaculatus \
--metadata ./Data/Halichoeres_trimaculatus/sample_metadata.csv \
--bam_files ./out/MK.F1_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.F1.Mix_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.F2_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.F3_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.F4_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M1_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M1.Mix_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M2_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M3_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M4_R1_trim.fq_Aligned.sortedByCoord.out.bam

```

&nbsp;

&nbsp;

&nbsp;

&nbsp;

## 20250414

- [./Tools/06_combined_exec_TSSr.cwl](./Tools/06_combined_exec_TSSr.cwl) での実行が完了
- 一応出力まで確認
- [./Tools/07_join_all_assignedClusters.cwl](./Tools/07_join_all_assignedClusters.cwl)での実行が完了
- 参考にRscriptの実行コマンドを掲載

&nbsp;

```bash
Rscript ./scripts/05-08_combined_exec_TSSr_modified.R \
--refSource ./Data/Halichoeres_trimaculatus/braker_correctID_v3.gtf \
--seedFile ./Data/Halichoeres_trimaculatus/BSgenome.Htrimaculatus.Htriv1.1-seed_2.txt \
--threads 16 \
--inputFilesType bamPairedEnd \
--organismName Halichoeres trimaculatus \
--metadata ./Data/Halichoeres_trimaculatus/sample_metadata.csv \
--upstream 1000 \
--downstream 0 \
--annotationType genes \
--bam_files ./out/MK.F1_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.F1.Mix_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.F2_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.F3_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.F4_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M1_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M1.Mix_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M2_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M3_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M4_R1_trim.fq_Aligned.sortedByCoord.out.bam
```

&nbsp;

```bash
cwltool --debug --outdir ./out/ ./Tools/06_combined_exec_TSSr.cwl ./config/Commandlinetool_config/06_tssr_config.yaml
cwltool --outdir ./out/ ./Tools/07_join_all_assignedClusters.cwl ./config/Commandlinetool_config/07_join_all_assignedClusters.yml 
```

&nbsp;

### TSSrの改良について

- 野津さんが現在TSSrの改良を進めている
- [RyoNozu/TSSr](https://github.com/RyoNozu/TSSr) 
- まだコンテナ化していないが、もしかしたら後日するかも

&nbsp;

### 今後の展望

- 明日サブワークフローを構築、できたらメインのワークフローに統合

&nbsp;

&nbsp;

&nbsp;

&nbsp;

## 20250415 memo by yonezawa

### TSSr のDocker image

- 1.0.2 からは野津さんが修正したTSSr (https://github.com/RyoNozu/TSSr) を実行するように変更
