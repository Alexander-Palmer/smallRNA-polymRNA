#!/bin/bash

#Removes RNA adapter (index) from different samples
#############################################################################

#Sample A-284-1_intestine_ALG1
cutadapt -a (TGGAATTCTCGGGTGCCAAGG \
         -o Cutadapt_out/A-284-1_intestine_ALG1_clip.fastq Raw_seq/A-284-1_intestine_ALG1.fastq\
         --cores=6
#Sample A-284-2_intestine_ALG1
cutadapt -a (TGGAATTCTCGGGTGCCAAGG \
         -o Cutadapt_out/A-284-2_intestine_ALG1_clip.fastq Raw_seq/A-284-2_intestine_ALG1.fastq\
         --cores=6

#Sample A-285-1_intestine_ALG2
cutadapt -a (TGGAATTCTCGGGTGCCAAGG \
         -o Cutadapt_out/A-285-1_intestine_ALG2_clip.fastq Raw_seq/A-285-1_intestine_ALG2.fastq\
         --cores=6
#Sample A-285-2_intestine_ALG2
cutadapt -a (TGGAATTCTCGGGTGCCAAGG \
         -o Cutadapt_out/A-285-2_intestine_ALG2_clip.fastq Raw_seq/A-285-2_intestine_ALG2.fastq\
         --cores=6


#Sample A-286-1_neuron_ALG1
cutadapt -a (TGGAATTCTCGGGTGCCAAGG \
         -o Cutadapt_out/A-286-1_neuron_ALG1_clip.fastq Raw_seq/A-286-1_neuron_ALG1.fastq\
         --cores=6
#Sample A-286-2_intestine_ALG1
cutadapt -a (TGGAATTCTCGGGTGCCAAGG \
         -o Cutadapt_out/A-286-2_neuron_ALG1_clip.fastq Raw_seq/A-286-2_neuron_ALG1.fastq\
         --cores=6

#Sample A-287-1_neuron_ALG2
cutadapt -a (TGGAATTCTCGGGTGCCAAGG \
         -o Cutadapt_out/A-287-1_neuron_ALG2_clip.fastq Raw_seq/A-287-1_neuron_ALG2.fastq\
         --cores=6
#Sample A-287-2_neuron_ALG2
cutadapt -a (TGGAATTCTCGGGTGCCAAGG \
         -o Cutadapt_out/A-287-2_neuron_ALG2_clip.fastq Raw_seq/A-287-2_neuron_ALG2.fastq\
         --cores=6


#Sample A-290-1_bwm_ALG1
cutadapt -a (TGGAATTCTCGGGTGCCAAGG \
         -o Cutadapt_out/A-290-1_bwm_ALG1_clip.fastq Raw_seq/A-290-1_bwm_ALG1.fastq\
         --cores=6
#Sample A-290-2_bwm_ALG1
cutadapt -a (TGGAATTCTCGGGTGCCAAGG \
         -o Cutadapt_out/A-290-2_bwm_ALG1_clip.fastq Raw_seq/A-290-2_bwm_ALG1.fastq\
         --cores=6

#Sample A-291-1_bwm_ALG2
cutadapt -a (TGGAATTCTCGGGTGCCAAGG \
         -o Cutadapt_out/A-291-1_bwm_ALG2_clip.fastq Raw_seq/A-291-1_bwm_ALG2.fastq\
         --cores=6
#Sample A-291-2_bwm_ALG2
cutadapt -a (TGGAATTCTCGGGTGCCAAGG \
         -o Cutadapt_out/A-291-2_bwm_ALG2_clip.fastq Raw_seq/A-291-2_bwm_ALG2.fastq\
         --cores=6

mv (*clip*) ../Cutadapt_out
cd ..
