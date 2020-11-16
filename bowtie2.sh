#!/bin/bash

#Align small RNA-seq data to C. elegans genome using bowtie2
############################################################

#Build bowtie2 index
bowtie2-build -f C.elegans.WBcel235.dna.toplevel.fa cel_index

#bowtie2 alignment
bowtie2 -x cel_index -U Intestine_ALG-1_rep1.fastq,Intestine_ALG-1_rep2.fastq,\
Neuron_ALG-1_rep1.fastq,Neuron_ALG-1_rep2.fastq,Muscle_ALG-1_rep1.fastq,\
Muscle_ALG-1_rep2.fastq,Intestine_ALG-2_rep1.fastq,Intestine_ALG-2_rep2.fastq,\
Neuron_ALG-2_rep1.fastq,Neuron_ALG-2_rep2.fastq,Muscle_ALG-2_rep1.fastq,\
Muscle_ALG-2_rep2.fastq -S seq
