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
