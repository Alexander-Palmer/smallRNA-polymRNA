#!/bin/bash

#Perform fastQC analysis on sequencing data before and after 
#adapter trimming
#####################################################################

fastqc S996-1_combined_R1.fastq S996-1_combined_R2.fastq \
      S996-2_combined_R1.fastq S996-2_combined_R2.fastq \
      S996-3_combined_R1.fastq S996-3_combined_R2.fastq \
      S997-1_combined_R1.fastq S997-1_combined_R2.fastq \
      S997-2_combined_R1.fastq S997-2_combined_R2.fastq \
      S997-3_combined_R1.fastq S997-3_combined_R2.fastq \
      S998-1_combined_R1.fastq S998-1_combined_R2.fastq \
      S998-2_combined_R1.fastq S998-2_combined_R2.fastq \
      S998-3_combined_R1.fastq S998-3_combined_R2.fastq \
      S999-1_combined_R1.fastq S999-1_combined_R2.fastq \
      S999-2_combined_R1.fastq S999-2_combined_R2.fastq \
      S999-3_combined_R1.fastq S999-3_combined_R2.fastq \
      S1000-1_combined_R1.fastq S1000-1_combined_R2.fastq \
      S1000-2_combined_R1.fastq S1000-2_combined_R2.fastq \
      S1000-3_combined_R1.fastq S1000-3_combined_R2.fastq \
      S1045-1_combined_R1.fastq S1045-1_combined_R2.fastq \
      S1045-2_combined_R1.fastq S1045-2_combined_R2.fastq \
      S1045-3_combined_R1.fastq S1045-3_combined_R2.fastq \
      S1046-1_combined_R1.fastq S1046-1_combined_R2.fastq \
      S1046-2_combined_R1.fastq S1046-2_combined_R2.fastq \
      S1046-3_combined_R1.fastq S1046-3_combined_R2.fastq \
      S1047-1_combined_R1.fastq S1047-1_combined_R2.fastq \
      S1047-2_combined_R1.fastq S1047-2_combined_R2.fastq \
      S1047-3_combined_R1.fastq S1047-3_combined_R2.fastq \
      -o ../FastQC_out -f fastq
