### DynamicEX introduction

---

DynamicEX is a single-cell transcriptome data analysis software that includes multiple analysis steps. The input file is the original sequencing Fastq file, and the output is the cell expression matrix for subsequent analysis, as well as a web report that displays various statistical indicators of the data.

### Install DynamicEX 1.0.0

----

#### Install dependencies

-----

The DynamicEX software package already contains most of the software required for analysis, and can be used after installing dependencies.

- Python (>= 3.6)
- Snakemake (== 5.4.5)
- Singularity (>= 3.7.2)

Python：

- Install Python from source code, visit [python](https://www.python.org/downloads/) for the download link.

- Using conda

```shell
conda create DynamicEX python==3.8
```

Snakemake:

- Using pip

```shell
pip install snakemake==5.4.5
```

- Using conda

```shell
conda install snakemake==5.4.5
```

Singularity:

- Installation must be done using root privilege. Please refer to the official installation guide at [singularity](https://github.com/sylabs/singularity/blob/main/INSTALL.md). RPM installation is recommended.

- Download sif files already built by DynamicEX. Download the compressed file from [here](https://bioservice.obs.cn-jssz1.ctyun.cn:443/hub%2Fsif.tar.gz?AccessKeyId=ZQMMROORLBRHRFUORXSZ&Expires=1694306336&Signature=7Lmyhd%2BaEm43DVEueos8M6aBdHY%3D).

```shell
cd DynamicEX
mkdir sif
tar -zxf hub_sif.tar.gz -C sif
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
  - --library: library version,default Decoder1.0
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





