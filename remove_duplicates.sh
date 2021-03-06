#!/bin/bash
#
#PBS -A qris-uq
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=10:00:00
#PBS -N Remove_duplicate_reads_6
#

#Remove gene duplicates generated during hisat2
###############################################

grep -vf gene_duplicates.txt S998-1.sam > S998-1_nodup.sam
grep -vf gene_duplicates.txt S998-2.sam > S998-2_nodup.sam
grep -vf gene_duplicates.txt S998-3.sam > S998-3_nodup.sam

grep -vf gene_duplicates.txt S1000-1.sam > S1000-1_nodup.sam
grep -vf gene_duplicates.txt S1000-2.sam > S1000-2_nodup.sam
grep -vf gene_duplicates.txt S1000-3.sam > S1000-3_nodup.sam
