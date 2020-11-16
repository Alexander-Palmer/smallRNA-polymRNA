#!/bin/bash
#
#PBS -A qris-uq
#PBS -l select=1:ncpus=8:mem=60GB
#PBS -l walltime=01:00:00
#PBS -N Sort_files_7_fast
#

#Sorts .bam files for featureCounts annotation
##############################################

samtools sort  S998-1-filtered.bam > S998-1-filtered-sorted.bam -@ 96 
samtools sort  S998-2-filtered.bam > S998-2-filtered-sorted.bam -@ 96 
samtools sort  S998-3-filtered.bam > S998-3-filtered-sorted.bam -@ 96 

samtools sort  S1000-1-filtered.bam > S1000-1-filtered-sorted.bam -@ 96 
samtools sort  S1000-2-filtered.bam > S1000-2-filtered-sorted.bam -@ 96 
samtools sort  S1000-3-filtered.bam > S1000-3-filtered-sorted.bam -@ 96 
