# Upstream analysis of RNAseq

---

- [Upstream analysis of RNAseq](#upstream-analysis-of-rnaseq)
  - [DOWNLOAD TOOLS](#download-tools)
    - [Setup Tools dir](#setup-tools-dir)
    - [SRA-Toolkit](#sra-toolkit)
    - [Download FASTQC](#download-fastqc)
    - [Download Trimmomatic](#download-trimmomatic)
    - [Download STAR](#download-star)
    - [Download subRead, featureCount](#download-subread-featurecount)
    - [Download RSeQC](#download-rseqc)
    - [Download Picard](#download-picard)
  - [SETUP WORKING DIR](#setup-working-dir)
    - [Results](#results)
    - [Download genome of S.cerevisiae](#download-genome-of-scerevisiae)
    - [Download gtf file](#download-gtf-file)
    - [Download BED file](#download-bed-file)
  - [PATH](#path)
  - [RNASEQ: RAW DATA PROCESSING \& ALIGNMENT](#rnaseq-raw-data-processing--alignment)
    - [Create a list of SRA Accession Number](#create-a-list-of-sra-accession-number)
    - [Then copy and paste these Acessions to your SraAccList.csv](#then-copy-and-paste-these-acessions-to-your-sraacclistcsv)
  - [COMMAND TO DOWNLOAD RAW DATA OF THE PROJECT](#command-to-download-raw-data-of-the-project)
  - [CHANGE THE FILE NAME](#change-the-file-name)
    - [Result](#result)
  - [RAW DATA PROCESSING: QUALITY CONTROL](#raw-data-processing-quality-control)
  - [RAW DATA PROCESSING: TRIMMING \& FILTERING](#raw-data-processing-trimming--filtering)
  - [ALIGNMENT WITH STAR](#alignment-with-star)
    - [Index the genome](#index-the-genome)
    - [Alignment](#alignment)
    - [Move all bam file to a location](#move-all-bam-file-to-a-location)
  - [POST-PROCESSING ALIGNMENT DATA: REMOVE DUPLICATE](#post-processing-alignment-data-remove-duplicate)
  - [POST-PROCESSING ALIGNMENT DATA: QUALITY CONTROL USING RSEQQC](#post-processing-alignment-data-quality-control-using-rseqqc)
    - [RSeqQC](#rseqqc)
  - [QUANTIFICATION](#quantification)
    - [featureCounts](#featurecounts)

## DOWNLOAD TOOLS

### Setup Tools dir

```bash
mkdir tools
```

### SRA-Toolkit

```bash
sudo apt update
sudo apt install sra-toolkit
fastq-dump --version
```

### Download FASTQC

```bash
## Install & unzip
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
## Write at the end of the .bashrc file. This command will let you export the path to FastQC, and
execute it everywhere
nano ~/.bashrc
export PATH='path/to/FastQC/':$PATH
## For example:
export PATH='/home/duydao/dnaseq_work/tools/FastQC':$PATH
source ~/.bashrc
## Try it by running
fastqc
```

### Download Trimmomatic

```bash
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
```

### Download STAR

```bash
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
tar -xzf 2.7.10b.tar.gz
cd STAR-2.7.10b

# Compile
cd source
make STAR

# Export path to STAR
nano ~/.bashrc
export PATH='path/to/tools/STAR-2.7.10b/bin/Linux_x86_64/':$PATH
```

### Download subRead, featureCount

```bash
git clone https://github.com/ShiLab-Bioinformatics/subread.git

cd subread/src
make -f Makefile.Linux

nano ~/.bashrc
export PATH='path/to/tools/subread/bin':$PATH
```

### Download RSeQC

Download at this link:
<https://sourceforge.net/projects/rseqc/files/RSeQC-5.0.1.tar.gz/download>

```bash
mv ~/Downloads/RSeQC-5.0.1.tar.gz tools/
gunzip RSeQC-5.0.1.tar.gz
tar -xvf RSeQC-5.0.1.tar.gz
```

### Download Picard

```bash
git clone https://github.com/broadinstitute/picard.git
cd picard/
./gradlew shadowJar
```

## SETUP WORKING DIR

```bash
mkdir sacCer2
mkdir -p sacCer2/sra/
mkdir -p sacCer2/ref/
mkdir -p sacCer2/ref/annotation/
mkdir -p sacCer2/ref/genome/
mkdir -p sacCer2/raw/
mkdir -p sacCer2/raw/qc_check
mkdir -p sacCer2/trim/
mkdir -p sacCer2/trim/qc_check
mkdir -p sacCer2/align
mkdir -p sacCer2/align/bam
mkdir -p sacCer2/counts/
```

### Results

```text
.
├── align
├── counts
├── raw
│   └── qc_check
├── ref
│   ├── annotation
│   └── genome
├── sra
└── trim
    └── qc_check
```

### Download genome of S.cerevisiae

```bash
cd $p_ref/genome/
wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz
```

### Download gtf file

```bash
cd $p_ref/annotation/
wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz
gunzip sacCer3.ensGene.gtf.gz
```

### Download BED file

Follow the guide at this link: <https://www.biostars.org/p/465274/>

Download the BED of sacCer3 at this link: <https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1648314330_jwKWTJMxFAfhUJN4aLNKfho0wfgN&clade=other&org=S.+cerevisiae&db=sacCer3&hgta_group=genes&hgta_track=ensGene&hgta_table=0&hgta_regionType=genome&position=chrIV%3A765%2C966-775%2C965&hgta_outputType=primaryTable&hgta_outFileName=>

## PATH

```bash
p_sra='your/path/to/sra'
p_raw='your/path/to/raw'
p_trim="your/path/to/trim"
trimmomatic_jar='your/path/to/trimmomatic-0.39.jar'
p_ref='your/path/to/ref'
p_align='your/path/to/align'
p_counts='your/path/to/counts'
p_bam='your/path/to/bam'
```

## RNASEQ: RAW DATA PROCESSING & ALIGNMENT

### Create a list of SRA Accession Number

```bash
cd $p_sra
nano SraAccList.csv
```

### Then copy and paste these Acessions to your SraAccList.csv

```text
SRR23867632
SRR23867633
SRR23867634
SRR23867635
SRR23867636
SRR23867637
```

## COMMAND TO DOWNLOAD RAW DATA OF THE PROJECT

```bash
for sra_acc in $(cat $p_sra/SraAccList.csv); do

 # Get the SRA
 prefetch -v $sra_acc --output-directory $p_sra
    echo "=== $sra_acc sra file downloaded ==="
 # Download fastq file
 fastq-dump \
    --outdir $p_raw/ \
    --split-files $p_sra/${sra_acc}/${sra_acc}.sra \
    && gzip -f $p_raw/*.fastq
    echo "=== Split $sra_acc file completed ==="

done
```

## CHANGE THE FILE NAME

```bash
cd $p_raw
mv SRR23867633_1.fastq.gz WT_C_1_R1.fastq.gz
mv SRR23867633_2.fastq.gz WT_C_1_R2.fastq.gz
mv SRR23867634_1.fastq.gz WT_E_1_R1.fastq.gz
mv SRR23867634_2.fastq.gz WT_E_1_R2.fastq.gz
mv SRR23867635_1.fastq.gz WT_C_2_R1.fastq.gz
mv SRR23867635_2.fastq.gz WT_C_2_R2.fastq.gz
mv SRR23867632_1.fastq.gz WT_E_2_R1.fastq.gz
mv SRR23867632_2.fastq.gz WT_E_2_R2.fastq.gz
mv SRR23867636_1.fastq.gz WT_C_3_R1.fastq.gz
mv SRR23867636_2.fastq.gz WT_C_3_R2.fastq.gz
mv SRR23867637_1.fastq.gz WT_E_3_R1.fastq.gz
mv SRR23867637_2.fastq.gz WT_E_3_R2.fastq.gz
```

### Result

```bash
tree -h raw/

├── [4.0K]  qc_check
├── [450M]  WT_C_1_R1.fastq.gz
├── [448M]  WT_C_1_R2.fastq.gz
├── [445M]  WT_C_2_R1.fastq.gz
├── [449M]  WT_C_2_R2.fastq.gz
├── [463M]  WT_C_3_R1.fastq.gz
├── [468M]  WT_C_3_R2.fastq.gz
├── [472M]  WT_E_1_R1.fastq.gz
├── [472M]  WT_E_1_R2.fastq.gz
├── [493M]  WT_E_2_R1.fastq.gz
├── [560M]  WT_E_2_R2.fastq.gz
├── [485M]  WT_E_3_R1.fastq.gz
└── [483M]  WT_E_3_R2.fastq.gz
```

## RAW DATA PROCESSING: QUALITY CONTROL

Use FastQC

```BASH
for file in $(ls $p_raw/*gz); do
    fastqc $file -o $p_raw/qc_check/
done
```

## RAW DATA PROCESSING: TRIMMING & FILTERING

Trim and remove bad quality sequence with Trimmomatic

```bash
## Loop over the files ending with ".gz"
cd $p_raw
for file in $(ls *_R1*.gz);do
    # Extract the base name without the extension
    base_name="${file%_*.*.*}"

    # Construct the input and output file names
    read_1=${base_name}_R1.fastq.gz
    read_2=${base_name}_R2.fastq.gz
    output_paired_1="${base_name}_R1_paired.fastq.gz"
    output_unpaired_1="${base_name}_R1_unpaired.fastq.gz"
    output_paired_2="${base_name}_R2_paired.fastq.gz"
    output_unpaired_2="${base_name}_R2_unpaired.fastq.gz"
    
    # Run Trimmomatic with the specified parameters
    java -jar $trimmomatic_jar PE \
        -threads 14 \
        -trimlog "$p_trim/trim_test.log" \
        "$p_raw/$read_1" \
        "$p_raw/$read_2" \
        "$p_trim/$output_paired_1" \
        "$p_trim/$output_unpaired_1" \
        "$p_trim/$output_paired_2" \
        "$p_trim/$output_unpaired_2" \
        HEADCROP:10 \
        CROP:65 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:30 \
        MINLEN:36
    
    # Remove unpaired file
    rm $p_trim/$output_unpaired_1 $p_trim/$output_unpaired_2

    # Check fastqc again
    fastqc $p_trim/$output_paired_1 -o $p_trim/qc_check
    fastqc $p_trim/$output_paired_2 -o $p_trim/qc_check
done
```

## ALIGNMENT WITH STAR

### Index the genome

```bash
# Using STAR to index the genome
STAR \
--runThreadN 12 \
--runMode genomeGenerate \
--genomeDir $p_ref/genome \
--genomeFastaFiles $p_ref/genome/sacCer3.fa
```

### Alignment

```bash
cd $p_trim 

# 
for file in $(ls *_R1*.gz);do
    # Extract the base name without the extension
    base_name="${file%_*_*.*.*}"

    # Construct the input and output file names
    read_1=${base_name}_R1_paired.fastq.gz
    read_2=${base_name}_R2_paired.fastq.gz

    # Align with STAR
    STAR \
    --runThreadN 12 \
    --runMode alignReads \
    --readFilesType Fastx \
    --readFilesCommand zcat \
    --genomeDir $p_ref/genome \
    --sjdbGTFfile $p_ref/annotation/sacCer3.ensGene.gtf \
    --sjdbOverhang 100 \
    --readFilesIn \
    $read_1 \
    $read_2 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMmode Full \
    --outFileNamePrefix $p_align/${base_name}_
done
```

### Move all bam file to a location

```bash
cd $p_align
mkdir bam
mv *.bam bam/
```

## POST-PROCESSING ALIGNMENT DATA: REMOVE DUPLICATE

```BASH
cd $p_align/bam/

# Create path to your picard tool
picard_jar='path/to/picard.jar'

# Using for loop to remove all duplicate in bam file
for file in $(ls *.bam);do
    base_name="${file%Aligned*}"
    java -jar $picard_jar MarkDuplicates \
    --INPUT ${base_name}Aligned.sortedByCoord.out.bam \
    --OUTPUT ${base_name}_rmdup.bam \
    --METRICS_FILE ${base_name}_rmdup.metrics2 \
    --REMOVE_DUPLICATES true
done
```

## POST-PROCESSING ALIGNMENT DATA: QUALITY CONTROL USING RSEQQC

Try to explore the quality of all of your BAM file using RSeqQC package

### RSeqQC

1. Index bam

    ```BASH
    samtools index WT_C_1_rmdup.bam
    ```

2. Genebody coverage

    ```BASH
    geneBody_coverage.py \
        -r $p_ref/annotation/sacCer3.ens.bed \
        -i $p_bam/WT_C_1_rmdup.bam  \
        -o gene_coverage
    ```

3. Junction annotation

    ```BASH
    junction_annotation.py \
        -r $p_ref/annotation/sacCer3.ens.bed \
        -i $p_bam/WT_C_1_rmdup.bam \
        -o WT_C_1
    ```

4. Junction saturation

    ```BASH
    junction_saturation.py \
        -r $p_ref/annotation/sacCer3.ens.bed \
        -i $p_bam/WT_C_1_rmdup.bam \
        -o WT_C_1
    ```

5. Read Distribution

    ```BASH
        read_distribution.py \
        -i $p_bam/WT_C_1_rmdup.bam \
        -r $p_ref/annotation/sacCer3.ens.bed
    ```

## QUANTIFICATION

### featureCounts

featureCounts counts how many reads map to genomic features, such as genes, exon, promoter and genomic bins.

```bash
featureCounts \
-p \
-a $p_annotation/sacCer3.ensGene.gtf \
-o $p_counts/all.featureCounts.txt \
$p_bam/WT_C_1_rmdup.bam \
$p_bam/WT_C_2_rmdup.bam \
$p_bam/WT_C_3_rmdup.bam \
$p_bam/WT_E_1_rmdup.bam \
$p_bam/WT_E_2_rmdup.bam \
$p_bam/WT_E_3_rmdup.bam
```
