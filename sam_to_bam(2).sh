#!/bin/bash
#
#PBS -A qris-uq
#PBS -l select=1:ncpus=16:mem=16GB
#PBS -l walltime=01:00:00
#PBS -N Sam_to_bam_8
#

#Convert hisat2 .sam output to .bam format
##########################################

samtools view -S -b S998-1_nodup.sam > S998-1.bam -@ 192 -m 16G
samtools view -S -b S998-2_nodup.sam > S998-2.bam -@ 192 -m 16G
samtools view -S -b S998-3_nodup.sam > S998-3.bam -@ 192 -m 16G

samtools view -S -b S1000-1_nodup.sam > S1000-1.bam -@ 192 -m 16G
samtools view -S -b S1000-2_nodup.sam > S1000-2.bam -@ 192 -m 16G
samtools view -S -b S1000-3_nodup.sam > S1000-3.bam -@ 192 -m 16G

mv *.bam ../bam_out
