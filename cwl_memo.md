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

&nbsp;

## 2025/02/09 memo

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