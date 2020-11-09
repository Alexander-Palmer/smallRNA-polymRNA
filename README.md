# smallRNA - polysome-seq
Scripts used to process and analyse small-RNA and polysome sequencing data.

## Dependencies
* fastQC (vx.x.x)
* cutadapt (v2.10)
* HISAT2 (v2.10)
* samtools (vx.x.x)
* FeatureCounts Subread (v2.0.0)
* edgeR (v3.11)
* biomaRt (v2.46.0)
 

## Pipeline
### smallRNA-seq
1. **Cutadapt**: La
2. **fastQC**: De
3. **bowtie2**: Da

### Polysome-seq

HISAT gene duplicates:
_C. elegans_ genome annotation file:
_C. elegans_ reference genome:


1. **fastQC.sh** - _initial assessment of read quality_
2. **cutadapt.sh** - _removal of sequencing adapters_
3. **fastQC.sh** - _assessment of read quality folowing adapter removal_
4. **hisat2.pbs** - _alignment of trimmed reads to the_ C. elegans _genome_
5.  **remove_duplicates.pbs** - _removal of gene duplicates generated during HISAT2 alignment_
6. **sam_to_bam.pbs** - _conversion of HISAT2 .sam output to .bam format_
7.  **filter_pairs.pbs** - _filtering of processed reads to include matching pairs only_
8.  **sort_bam.pbs** - _sorting of .bam files for annotation_
9.  **featureCounts_annotation.R** - _annotation and quantification of aligned reads_
10. **edgeR_polysome_analysis.R** - _normalisation and analysis of mRNAs, including multiple comparisons_



## Citation
