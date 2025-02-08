# CWLを書く手順

## 1. zatsu-cwl-generatorでざっくりCWLファイルを出力する

```bash
zatsu-cwl-generator "seqkit ./Data/Halichoeres_trimaculatus/Halichoeres_trimaculatus-hifiasm-3ddna-v1.1.edit.fna.gz -i -O ./seqs_by_scaffold" > ./Tools/01_split_genome_seqs.cwl
```

## 2. 手動で修正，チェックする

- CWLの記法のバージョンを1.2にする
- 以下のコマンドでvalidateして書き方に変なところがないか確認する

```bash
cwltool --validate ./Tools/01_split_genome_seqs.cwl 
```

## 3. チェックが通ったら，ワークフローエンジンのcwltoolで実行する

```bash
cwltool --debug --cachedir ./cwl_cache/ --outdir ./Data/ ./Tools/01_split_genome_seqs.cwl
```

&nbsp;

&nbsp;

&nbsp;

## 2025/02/08 memo

- fastpのプロセスで複数ファイルの処理を行うCWLファイルを新たに作成
- サブワークフローを作成して，複数ファイルの処理を行う
- YAMLファイルにファイルのパスを記載して，サブワークフローを実行


```bash
cwltool --debug --outdir ./out/ ./workflow/01_trimming_fastq_subworkflow.cwl ./config/01_trimming_fastq_files.yml
```



