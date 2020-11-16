#!/bin/bash

#Convert bowtie2 .sam output to .fa format for miRDeep2 analysis
################################################################

for file in *.sam
do
  samtools bam2fq $file | seqtk seq -A > $file.fa
