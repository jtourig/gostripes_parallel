# GOSTRIPEs v2.0.0

This is a reduced and efficient implementation of GOSTRIPEs for the processing of STRIPE-seq data.
Only paired-end reads are currently supported.

## Installation

Clone the repository to a clean working directory that will be used for the analysis.

```
git clone https://github.com/rpolicastro/gostripes_parallel.git
```

It is recommended to use conda to install the required software.
Conda installs all software and dependencies to distinct software environments.
These environments act independently of system software,
and help avoid conflicts in dependencies.
Environments are also self contained,
so removing them if desired is easy.
Installation instructions for miniconda3 can be found [here](https://docs.conda.io/en/latest/miniconda.html).

After miniconda3 is installed, create a GOSTRIPEs environment.

```
conda create -n GOSTRIPEs -y -c conda-forge -c bioconda \
fastqc multiqc csvtk star samtools cutadapt parallel seqtk umi_tools
```

The GOSTRIPEs software environment can now be activated for analysis.

```
conda activate GOSTRIPEs
```

## Useage

GOSTRIPEs requires a genome assembly in FASTA format,
a genome annotation in GTF format,
and paired-end FASTQ files from STRIPE-seq.

For this example we'll use paired end S. cerevisiae (Sc) S288C STRIPE-seq data.

```
if [ ! -d seqs ]; then mkdir seqs; fi

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR107/019/SRR10759419/SRR10759419_1.fastq.gz | \
  gunzip > seqs/SRR10759419_1.fastq

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR107/019/SRR10759419/SRR10759419_2.fastq.gz | \
  gunzip > seqs/SRR10759419_2.fastq
```

We also need a comma separated (CSV) sample sheet with columns 'name', 'fastq_1', and 'fastq_2'
corresponding to the sample name, R1 FASTQ file, and R2 FASTQ file.

```
cat \
  <(echo name,fastq_1,fastq_2) \
  <(echo S288C,seqs/SRR10759419_1.fastq,seqs/SRR10759419_2.fastq) \
  > sample_sheet.csv
```

Finally, an Sc FASTA genome assembly and GTF genome annotation file will be downloaded from ENSEMBL.

```
if [ ! -d genome ]; then mkdir genome; fi

curl ftp://ftp.ensembl.org/pub/release-102/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.102.gtf.gz | \
  gunzip > genome/Saccharomyces_cerevisiae.R64-1-1.102.gtf

curl ftp://ftp.ensembl.org/pub/release-102/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz | \
  gunzip > genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa
```

Now that we have all of the required files, we can run the workflow.

```
bash gostripes.sh \
  -s sample_sheet.csv \
  -a genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa \
  -g genome/Saccharomyces_cerevisiae.R64-1-1.102.gtf \
  -t 4 \
  --paired \
  --genomeSAindexNbases 10
```
