#!/bin/bash

#Perform fastQC analysis on sequencing data prior to and after 
#cutadapt(2) adapter trimming
##############################################################

#Before cutadapt(2) adapter trimming
fastqc S998-1_combined_R1.fastq S998-1_combined_R2.fastq \
      S998-2_combined_R1.fastq S998-2_combined_R2.fastq \
      S998-3_combined_R1.fastq S998-3_combined_R2.fastq \
      S1000-1_combined_R1.fastq S1000-1_combined_R2.fastq \
      S1000-2_combined_R1.fastq S1000-2_combined_R2.fastq \
      S1000-3_combined_R1.fastq S1000-3_combined_R2.fastq \
      -o ../FastQC_out -f fastq

#After cutadapt(2) adapter trimming
fastqc S998-1_R1.fastq S998-1_R2.fastq \
      S998-2_R1.fastq S998-2_R2.fastq \
      S998-3_R1.fastq S998-3_R2.fastq \
      S1000-1_R1.fastq S1000-1_R2.fastq \
      S1000-2_R1.fastq S1000-2_R2.fastq \
      S1000-3_R1.fastq S1000-3_R2.fastq \
      -o ../FastQC_out -f fastq
