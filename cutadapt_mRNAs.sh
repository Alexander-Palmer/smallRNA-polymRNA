#!/bin/bash

#Removes RNA adapter (index) from different samples
#############################################################################

#Sample S998-1-R1
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG\
         -o Cutadapt_out/S998-1-R1.fastq Raw_seq/S998-1_combined_R1.fastq\
         --cores=6
#Sample S998-1-R2
cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT\
         -o Cutadapt_out/S998-1-R2.fastq Raw_seq/S998-1_combined_R2.fastq\
         --cores=6
#Sample S998-2-R1
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG\
         -o Cutadapt_out/S998-2-R1.fastq Raw_seq/S998-2_combined_R1.fastq\
         --cores=6
#Sample S998-2-R2
cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT\
         -o Cutadapt_out/S998-2-R2.fastq Raw_seq/S998-2_combined_R2.fastq\
         --cores=6
#Sample S998-3-R1
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG\
         -o Cutadapt_out/S998-3-R1.fastq Raw_seq/S998-3_combined_R1.fastq\
         --cores=6
#Sample S998-3-R2
cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT\
         -o Cutadapt_out/S998-3-R2.fastq Raw_seq/S998-3_combined_R2.fastq\
         --cores=6


#Sample S1000-1-R1
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG\
         -o Cutadapt_out/S1000-1-R1.fastq Raw_seq/S1000-1_combined_R1.fastq\
         --cores=6
#Sample S1000-1-R2
cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT\
         -o Cutadapt_out/S1000-1-R2.fastq Raw_seq/S1000-1_combined_R2.fastq\
         --cores=6
#Sample S1000-2-R1
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG\
         -o Cutadapt_out/S1000-2-R1.fastq Raw_seq/S1000-2_combined_R1.fastq\
         --cores=6
#Sample S1000-2-R2
cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT\
         -o Cutadapt_out/S1000-2-R2.fastq Raw_seq/S1000-2_combined_R2.fastq\
         --cores=6
#Sample S1000-3-R1
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG\
         -o Cutadapt_out/S1000-3-R1.fastq Raw_seq/S1000-3_combined_R1.fastq\
         --cores=6
#Sample S1000-3-R2
cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT\
         -o Cutadapt_out/S1000-3-R2.fastq Raw_seq/S1000-3_combined_R2.fastq\
         --cores=6

mv !(*combined*) ../Cutadapt_out
cd ..
