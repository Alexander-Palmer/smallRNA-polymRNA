# smallRNA - polysome-seq
Scripts used to process and analyse small-RNA and polysome sequencing data.

## Dependencies
* bowtie2 (v2.3.4.3)
* samtools (v1.11)
* miRDeep2 (v2.0.1.2)
* fastQC (v0.11.9)
* cutadapt (v2.10)
* HISAT2 (v2.10)
* FeatureCounts (Subread, v2.0.0)
* edgeR (v3.11)
* biomaRt (v2.46.0)
 

## Pipeline
### smallRNA-seq
1. **fastQC(1).sh**: _initial assessment of read quality_
2. **cutadapt(1).sh**: _removal of sequencing adapters_
3. **fastQC(1).sh**: _assessment of read quality folowing adapter removal_
4. **bowtie2.sh**: _index and alignment to the_ C. elegans _genome_
5. **sam_to_bam.sh**: _conversion of bowtie2 .sam output to .bam format_
5. **miRDeep2.sh**: _pre-processing and identification novel miRNAs_ 

### Polysome-seq

HISAT gene duplicates: **gene_duplicates.txt**\
_C. elegans_ genome annotation: **C.elegans.WBCel235.96.gtf.gz**\
_C. elegans_ reference genome: http://ensembl.org/info/data/ftp/index.html\


1. **fastQC(2).sh** - _initial assessment of read quality_
2. **cutadapt(2).sh** - _removal of sequencing adapters_
3. **fastQC(2).sh** - _assessment of read quality folowing adapter removal_
4. **hisat2.pbs** - _alignment of trimmed reads to the_ C. elegans _genome_
5.  **remove_duplicates.pbs** - _removal of gene duplicates generated during HISAT2 alignment_
6. **sam_to_bam.pbs** - _conversion of HISAT2 .sam output to .bam format_
7.  **filter_pairs.pbs** - _filtering of processed reads to include matching pairs only_
8.  **sort_bam.pbs** - _sorting of .bam files for annotation_
9.  **featureCounts_annotation.R** - _annotation and quantification of aligned reads_
10. **edgeR_polysome_analysis.R** - _normalisation and analysis of mRNAs, including multiple comparisons_



## Citation
