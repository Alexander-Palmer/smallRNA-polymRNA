#!/bin/bash

#Conduct miRDeep2 pre-processing and analysis on small 
RNA-seq data to identify and quantify novel miRNAs
######################################################

#Pre-procesing sequencing data
mapper.pl *seq.fa -c -m -l 18 -p C.elegans.WBcel235.dna.toplevel.fa\
-s reads_collapsed.fa -t reads_vs_ref_genome.arf -v

#Utilising main miRDeep2 algorithm
miRDeep2.pl -b 3 reads_collapsed.fa C.elegans.WBcel235.dna.toplevel.fa\
reads_vs_ref_genome.arf miRBase_WBcel235_mature.fa -t cel 2>report.log
