#!/bin/bash

#Align small RNA-seq data to C. elegans genome using bowtie2
############################################################

#Build bowtie2 index
bowtie2-build -f C.elegans.WBcel235.dna.toplevel.fa cel_index

#bowtie2 alignment
bowtie2 -x cel_index -U \
A-284-1_intestine_ALG1_clip.fastq,A-284-2_intestine_ALG1_clip.fastq,\
A-285-1_intestine_ALG2_clip.fastq,A-285-2_intestine_ALG2_clip.fastq,\
A-286-1_neuron_ALG1_clip.fastq,A-286-2_neuron_ALG1_clip.fastq,\
A-287-1_neuron_ALG2_clip.fastq,A-287-2_neuron_ALG2_clip.fastq,\
A-290-1_bwm_ALG1_clip.fastq,A-290-2_bwm_ALG1_clip.fastq,\
A-291-1_bwm_ALG2_clip.fastq,A-291-2_bwm_ALG2_clip.fastq,\
N2-1_clip.fastq,N2-2_clip.fastq\
-S mapped_seq
