#!/bin/bash

#Convert .sam bowtie2 output to .bam format for further processing and analysis
###############################################################################

#Intestine ALG-1 & ALG-2
samtools view -b A-284-1_intestine_ALG1_clip.sam > A-284-1_intestine_ALG1_clip.bam
samtools view -b A-284-2_intestine_ALG1_clip.sam > A-284-2_intestine_ALG1_clip.bam
samtools view -b A-285-1_intestine_ALG2_clip.sam > A-285-1_intestine_ALG2_clip.bam
samtools view -b A-285-2_intestine_ALG2_clip.sam > A-285-2_intestine_ALG2_clip.bam

#Neuron ALG-1 & ALG-2
samtools view -b A-286-1_neuron_ALG1_clip.sam > A-286-1_neuron_ALG1_clip.bam
samtools view -b A-286-2_neuron_ALG1_clip.sam > A-286-2_neuron_ALG1_clip.bam
samtools view -b A-287-1_neuron_ALG2_clip.sam > A-287-1_neuron_ALG2_clip.bam
samtools view -b A-287-2_neuron_ALG2_clip.sam > A-287-2_neuron_ALG2_clip.bam

#Muscle ALG-1 & ALG-2
samtools view -b A-290-1_bwm_ALG1_clip.sam > A-290-1_bwm_ALG1_clip.bam
samtools view -b A-290-2_bwm_ALG1_clip.sam > A-290-2_bwm_ALG1_clip.bam
samtools view -b A-291-1_bwm_ALG2_clip.sam > A-291-1_bwm_ALG2_clip.bam
samtools view -b A-291-2_bwm_ALG2_clip.sam > A-291-2_bwm_ALG2_clip.bam
