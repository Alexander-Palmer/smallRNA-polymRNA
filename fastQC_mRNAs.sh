#!/bin/bash

#Perform fastQC analysis on sequencing data prior to and after 
#adapter trimming
##############################################################

#Before cutadapt adapter trimming
fastqc S998-1_combined_R1.fastq S998-1_combined_R2.fastq \
      S998-2_combined_R1.fastq S998-2_combined_R2.fastq \
      S998-3_combined_R1.fastq S998-3_combined_R2.fastq \
      S1000-1_combined_R1.fastq S1000-1_combined_R2.fastq \
      S1000-2_combined_R1.fastq S1000-2_combined_R2.fastq \
      S1000-3_combined_R1.fastq S1000-3_combined_R2.fastq \
      -o ../FastQC_out -f fastq

fastqc S998-1_R1.fastq S998-1_combined_R2.fastq \
      S998-2_combined_R1.fastq S998-2_combined_R2.fastq \
      S998-3_combined_R1.fastq S998-3_combined_R2.fastq \
      S1000-1_combined_R1.fastq S1000-1_combined_R2.fastq \
      S1000-2_combined_R1.fastq S1000-2_combined_R2.fastq \
      S1000-3_combined_R1.fastq S1000-3_combined_R2.fastq \
      -o ../FastQC_out -f fastq
