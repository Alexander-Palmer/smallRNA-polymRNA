#!/bin/bash
#
#PBS -A qris-uq
#PBS -l select=1:ncpus=16:mem=24GB
#PBS -l walltime=01:00:00
#PBS -N hisat_alignment
#

#Align trimmed reads to reference genome (hisat index)
######################################################

./hisat2 -x ../HISAT_index/hisat_index_2 -1 ../FASTQ_files/S998-1-R1.fastq -2 ../FASTQ_files/S998-1-R2.fastq -S /30days/uqapalm2/S998-1.sam -p 16
./hisat2 -x ../HISAT_index/hisat_index_2 -1 ../FASTQ_files/S998-2-R1.fastq -2 ../FASTQ_files/S998-2-R2.fastq -S /30days/uqapalm2/S998-2.sam -p 16
./hisat2 -x ../HISAT_index/hisat_index_2 -1 ../FASTQ_files/S998-3-R1.fastq -2 ../FASTQ_files/S998-3-R2.fastq -S /30days/uqapalm2/S998-3.sam -p 16

./hisat2 -x ../HISAT_index/hisat_index_2 -1 ../FASTQ_files/S1000-1-R1.fastq -2 ../FASTQ_files/S1000-1-R2.fastq -S /30days/uqapalm2/S1000-1.sam -p 16
./hisat2 -x ../HISAT_index/hisat_index_2 -1 ../FASTQ_files/S1000-2-R1.fastq -2 ../FASTQ_files/S1000-2-R2.fastq -S /30days/uqapalm2/S1000-2.sam -p 16
./hisat2 -x ../HISAT_index/hisat_index_2 -1 ../FASTQ_files/S1000-3-R1.fastq -2 ../FASTQ_files/S1000-3-R2.fastq -S /30days/uqapalm2/S1000-3.sam -p 16

mv *.sam ../hisat2_out
