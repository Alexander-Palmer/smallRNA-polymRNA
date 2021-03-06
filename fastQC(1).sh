#!/bin/bash

#Perform fastQC analysis on small RNA-seq data prior to and after 
#cutadapt(1) adapter trimming
#################################################################

#Before cutadapt(1) adapter clipping
fastqc A-284-1_intestine_ALG1.fastq A-284-2_intestine_ALG1.fastq \
       A-285-1_intestine_ALG2.fastq A-285-2_intestine_ALG2.fastq \
       A-286-1_neuron_ALG1.fastq A-286-2_neuron_ALG1.fastq \
       A-287-1_neuron_ALG2.fastq A-287-2_neuron_ALG2.fastq \
       A-290-1_bwm_ALG1.fastq A-290-2_bwm_ALG1.fastq \
       A-291-1_bwm_ALG2.fastq A-291-2_bwm_ALG2.fastq \
       N2-1.fastq N2-2.fastq \
       -o ../FastQC_out -f fastq

#After cutadapt(1) adapter clipping
fastqc A-284-1_intestine_ALG1_clip.fastq A-284-2_intestine_ALG1_clip.fastq \
       A-285-1_intestine_ALG2_clip.fastq A-285-2_intestine_ALG2_clip.fastq \
       A-286-1_neuron_ALG1_clip.fastq A-286-2_neuron_ALG1_clip.fastq \
       A-287-1_neuron_ALG2_clip.fastq A-287-2_neuron_ALG2_clip.fastq \
       A-290-1_bwm_ALG1_clip.fastq A-290-2_bwm_ALG1_clip.fastq \
       A-291-1_bwm_ALG2_clip.fastq A-291-2_bwm_ALG2_clip.fastq \
       N2-1_clip.fastq N2-2_clip.fastq \
       -o ../FastQC_out -f fastq
