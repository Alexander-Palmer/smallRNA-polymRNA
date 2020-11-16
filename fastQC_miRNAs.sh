#!/bin/bash

#Perform fastQC analysis on small RNA-seq data before and after 
#adapter trimming
###############################################################

fastqc A-284-1_intestine_ALG1.fastq A-284-2_intestine_ALG1.fastq \
       A-285-1_intestine_ALG2.fastq A-285-2_intestine_ALG2.fastq \
       A-286-1_neuron_ALG1.fastq A-286-2_neuron_ALG1.fastq \
       A-287-1_neuron_ALG2.fastq A-287-2_neuron_ALG2.fastq \
       A-290-1_bwm_ALG1.fastq A-290-2_bwm_ALG1.fastq \
       A-291-1_bwm_ALG2.fastq A-291-2_bwm_ALG2.fastq \
       -o ../FastQC_out -f fastq
