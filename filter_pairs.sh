#!/bin/bash
#
#PBS -A qris-uq
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=10:00:00
#PBS -N Filter_pairs
#

#Filter processed reads by those with matching pairs only
#########################################################

samtools view -bf 1 S998-1.bam > S998-1-filtered.bam -@ 192 -m 8G
samtools view -bf 1 S998-2.bam > S998-2-filtered.bam -@ 192 -m 8G
samtools view -bf 1 S998-3.bam > S998-3-filtered.bam -@ 192 -m 8G

samtools view -bf 1 S1000-1.bam > S1000-1-filtered.bam -@ 192 -m 8G
samtools view -bf 1 S1000-2.bam > S1000-2-filtered.bam -@ 192 -m 8G
samtools view -bf 1 S1000-3.bam > S1000-3-filtered.bam -@ 192 -m 8G
