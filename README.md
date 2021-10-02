# Book_Mastering-RNA-Seq-Data-Analysis-in-Infectious-Diseases-Practices
This book covers the skills necessary to analyses RNA seq data, such as preprocessing under a Linux shell environment, quality control assessment, and differential expression analysis using R/Bioconductor packages under RStudio

<img src="https://github.com/nzungizelab/Book_Mastering-RNA-Seq-Data-Analysis-in-Infectious-Diseases-Practices/blob/NzungizeL/images/Mastering%20RNA%20Seq%20Data.jpg" 
     width="400" height="500">

This repository contains the entire book : Mastering RNA Seq Data Analysis in Infectious Diseases Practices (Malaria and COVID-19)

# Table of Contents
* Chapter One: Case of Malaria in Africa 2014-2015
  * Introduction
  * Data reproducibility
  * Exercise design
  * Raw data
  * Reference genome annotation files
  * Quality Control
  * Filtering raw data (fastq files)
  * Alignment of paired end reads to the genome
  * Visualizing Mapped Reads with IGV
  * Quantification (counting reads)
  * Differential expression analysis using edgeR
  * Session Info
  * Reference
  * Resources used

* Chapter Two: Case of COVID-19 in 2020
  * Introduction
  * Data reproducibility
  * Kalisto-Sleuth Workflow
  * Prepare your workplace
  * Download and rename raw data
  * Reference transcripts for Homo sapiens
  * Trimming and quality control
  * Transcripts quantification with Kallisto
  * Quantification using Kallisto
  * Differential transcript expression analysis
  * Visualizations of the expression data using PCA
  * Creation of data models
  * Generate the results table for analysis
  * Plot the normalized bootstraps across all samples
  * Interactive graphical visualization
  * Session Info
  * Reference
  * Resources used


## Introduction
In the practical exercise, we will identify the expression levels of Anopheles funestus in the presence of insecticides (Pyrethroid). Sequencing was done with the Illumina platform which produces RNA reads. The goal of this hands-on session is to perform some basic tasks in the downstream analysis of RNA-seq data. We will check the quality of raw data using fastqc, index the reference with HISat2 then aligned it to the reference genome (Anopheles funestus, AFunF3) using HISat2. We will visualize the bam file with IGV, reads quantification, and generate the count matrix of the aligned reads using the featureCounts then perform differential expression analysis using edgeR.

## Data reproducibility
In this practical RNA seq analysis guide, we will use RNA seq data from Anopheles funestus
which is a major malaria vector in Africa.

* Data source: Original data from [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/PRJEB24351?show=reads) and accessible through accession number
`PRJEB24351`.
* Paper reference: Weedall et al., 2019 " [PMID: 30894503](https://pubmed.ncbi.nlm.nih.gov/30894503/).
* Samples: RNA seq data of mosquitoes collected from Cameroun and Ghana exposed to the insecticide called pyrethroid. We have three samples, and two conditions (two
populations of mosquitoes exposed to the insecticide and one population of mosquitoes no exposed used as control).
* Question: which genes are differentially expressed in insecticide resistance mosquitoes compared to the unexposed mosquitoes.

## Exercise design

* The samples used in the practice section for RNA seq analysis are summarized in the table below:

File name | Rename sample | Condition 
------------ | -------------| -------------
ERR2240088_1.fastq.gz, ERR2240088_2.fastq.gz | FANG-4_R1.fastq.gz, FANG-4_R2.fastq.gz | Control
ERR2240089_1.fastq.gz, ERR2240089_2.fastq.gz | GHA-PER-1_R1.fastq.gz, GHA-PER-1_R2.fastq.gz | Treated
ERR2240092_1.fastq.gz, ERR2240092_2.fastq.gz | CMR-PER-1_R1.fastq.gz, CMR-PER-1_R2.fastq.gz | Treated

* Notes :
1. First create the main directory "RNAseq_Weedall" and subdirectories for each step of RNA seq analysis.
2. All tools used in this practical exercise are open sources and freely available.
3. When you see a "#" starts or ends the line in the code, it's a comment/explain the purpose of the command line.
4. Navigate to your subdirectories using "cd" followed by the name of the directory you want.

* Tools used:
1. Quality control: `fastqc`
2. Alignment: `hisat2`
3. Quantification: `featureCounts`
4. Visualizing Mapped Reads: `IGV`
5. Differential expression analysis: `edgeR`

- Create a working main directory (Open terminal in Linux)
```
cd ~/home/lambert/
mkdir RNAseq_Weedall
```
- Create a subdirectory for raw data
```
mkdir 1.Raw_data
```
- Create other subdirectories 1 to 6 necessary for RNA-seq analysis workflow
```
1.Raw_data # Contain the downloaded raw data (fastq.gz)
2.Ref_genome # Contain the reference genome and genome annotation file (.GTF/.GFF)
3.Quality_control # Quality control of raw data.
4.Trimmed # Trimmed reads.
5.Alignment # Mapped reads (Bam files)
6.Quantification # Summarized gene counts across all samples.
7.DE # Differential expression analysis.
```

* Exercise design
In total, we have six pairs of fastq files, four pairs belonging to the treated (CMR-1 and GHA-1) condition, and two pairs belonging to the control (FANG-4).
  * CMR: Cameroon
  * GHA: Ghana
  * FANG: is a fully insecticide susceptible strain derived from Angola.

* Question 1: How to move to the directory you just created?
  * Open terminal in Linux
```
cd ~/RNAseq_Weedall/1.Raw_data
```
* Question 2: How can you check your current working directory?
```
pwd
```
* Output: (make sure we're in the right spot!)
```
/home/lambert/RNAseq_Weedall
```


* Note: you may have a different path to the RNAseq_Weedall directory on your own machine.

* Acquisition of Raw Data (downloads)
  * Time required: depend on the internet connection and the power of the computer.
  * Open the Terminal. First, go to the folder, where the data should be stored.
  * Once you’re in the correct folder, make a 1.Raw_data directory if not already made as shown above (preliminary exercise).
```
 cd ~/home/lambert/RNAseq_Weedall/
mkdir 1.Raw_data
cd ~/RNAseq_Weedall/1.Raw_data
```
* Link to downloading the raw data
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/008/ERR2240088/ERR2240088_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/008/ERR2240088/ERR2240088_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/009/ERR2240089/ERR2240089_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/009/ERR2240089/ERR2240089_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/002/ERR2240092/ERR2240092_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/002/ERR2240092/ERR2240092_2.fastq.gz
```

* Notes :
1. If you want all datasets to visit, click [here](https://www.ebi.ac.uk/ena/browser/home) and add PRJEB24351 in the search tool then select any file you want to download under the column of FASTQ FTP. Then rightclick on the file and copy the link address then paste in your terminal "wget copied link address".
2. Each data should be downloaded by its own FTP link.
3. Remember to keep raw data files and results in separate directories.
4. You will find samples (2 replicate in each of 7 samples), download the samples corresponding to 'FANG-4' as control, 'GHA-PER-1' and 'CMR-PER-1'.
5. Rename the files into simple names according to "Submitted FTP" e.g: CMR-PER-1_R1.fastq.gz, CMR-PER-1_R2.fastq.gz.
6. CMR: sample collected in Cameroon and sample collected in Ghana (GHA).
7. PER: Permethrin (insecticide belong to the class of pyrethroid)
8. CMR-PER-1: mosquitoes alive after 1 hour of PERmethrin exposure collected in Cameroon.

* Downloading reference genome
Visit the website of [vectorbase](https://vectorbase.org/vectorbase/app) in `Data` tab select `Download data files` and click Current_Release/ then in the list of species click `AfunestusFUMOZ/` directory next click `fasta/directory` click `Data` then download VectorBase-53_AfunestusFUMOZ_Genome.fasta

All commands that are given in this practical tutorial should be run within the main directory called `RNAseq_Weedall`.
```
~/home/lambert/RNAseq_Weedall/
```

* Raw data
Check that the raw data (paired-end reads) directory contains the above-mentioned files by typing:
```
cd ~/RNAseq_Weedall/1.Raw_data
ls -lh
```
 * Output: 
 
<img src="https://github.com/nzungizelab/Book_Mastering-RNA-Seq-Data-Analysis-in-Infectious-Diseases-Practices/blob/NzungizeL/images/1.png" 
     width="500" height="200">
* Alternative option :
```
ls -l *.fastq.gz
```
* Check the first 5 lines of reads:
```
zcat CMR-PER-1_R1.fastq.gz | head -n 5
```
* Alternative option :
```
head -5 CMR-PER-1_R1.fastq
```
* Check the last 4 lines of reads:
```
zcat CMR-PER-1_R1.fastq.gz | tail -n 4
```

* Output: 
 
<img src="https://github.com/nzungizelab/Book_Mastering-RNA-Seq-Data-Analysis-in-Infectious-Diseases-Practices/blob/NzungizeL/images/2.png" 
     width="500" height="200">
     
* We can view the file using the command called "less":
```
gzip -dc FANG-4_R1.fastq.gz | less
```
* Check and view the head of the file
```
zcat *.fastq.gz | head
```
