### DynamicEX introduction

---

DynamicEX is a single-cell transcriptome data analysis software that includes multiple analysis steps. The input file is the original sequencing Fastq file, and the output is the cell expression matrix for subsequent analysis, as well as a web report that displays various statistical indicators of the data.

### Install DynamicEX 1.0.0

----

Dowload DynamicEX from [here]()

```shell
mkdir DynamicEX

tar -zxf DynamicEX.tar.gz -C DynamicEX

source DynamicEX/bin/activate # This command needs to be executed every time DynamicEX is used

conda-unpack
```

### Manual

---

#### params

---

- mkref
  - --forceall: Force the execution of the selected (or the first) rule and all rules it is dependent on regardless of already created output.
  - --dryrun : Do not execute anything, and display what would be done. If you have a very large workflow, use --dryrun --quiet to just print a summary of the DAG of jobs.
  - --cores: set max cores the pipeline may request at one time.;default 16
  - --genomeName: output dir, used to build genome
  - --fasta: path to FASTA file containing your genome reference
  - --gtf: Path to genes GTF file containing annotated genes for your genome reference
  - --threads: Number of threads used during STAR genome index generation. Defaults to 8
- count
  - --sample: Sample name of each batch,required
  - --id: Final sample name,required
  - --inputdir: raw data path,required
  - --gtf: genome annotation file,required
  - --transcriptome: Path of folder containing transcriptome reference,required
  - --outputdir: output file name,default ./
  - --expect-cells: expect-cells,default 3000
  - --CBstart:  barcode start,default 26
  - --CBlen: barcode len,default 12
  - --UMIstart: umi start,default 38
  - --UMIlen: umi len,default 8
  - --cellcalling: method of cellcalling,default EmptyDrops_CR
  - --include_introns: Include intronic reads in count
  - --library: library version,default DECODER_1.0
  - --whitelist: use barcode whitelist

#### Quick start

---

- mkref

```shell
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.99.gtf.gz

DynamicEX mkref \
 --genome_name Homo_sapiens_GRCh38 \
 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
 --gtf Homo_sapiens.GRCh38.99.gtf
```

- count

```shell
#The input folder must contain fastq files, with the same file format as sample_S1_L001_R1_001.fastq.gz
DynamicEX count --sample sample --id S1 --inputdir rawdata/ --gtf Homo_sapiens.GRCh38.99.gtf --transcriptome Homo_sapiens_GRCh38 --outputdir result
```

The detailed documentation of the web_summary.html in the results can be found at [summary](https://github.com/DynamicBiosystems/DynamicEX/blob/main/doc/web_summary.md)

### About

---

DynamicEX is developed by dynamic-biosystems Co., Ltd. The official website is [www.dynamic-biosystems](http://www.dynamic-biosystems.com/).





