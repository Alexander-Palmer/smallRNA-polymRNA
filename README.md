# smallRNA - polysome-seq
Scripts used to process and analyse small-RNA and polysome sequencing data. Annotated using WBcel235.

## Dependencies
* fastQC (v0.11.9)
* cutadapt (v2.10)
* bowtie2 (v2.3.4.3)
* samtools (v1.11)
* miRDeep2 (v2.0.1.2)
* FeatureCounts (Subread, v2.0.0)
* edgeR (v3.11)
* biomaRt (v2.46.0)
* HISAT2 (v2.10)

## Pipeline
### smallRNA-seq

_C. elegans_ genome annotation: **C.elegans.WBCel235.96.gtf.gz** \
_C. elegans_ reference genome: http://ensembl.org/info/data/ftp/index.html \
Mature _C. elegans_ miRNAs: http://www.mirbase.org/cgi-bin/mirna_summary.pl?org=cel

1. **fastQC(1).sh**: _initial assessment of read quality_
2. **cutadapt(1).sh**: _removal of sequencing adapters_
3. **fastQC(1).sh**: _assessment of read quality folowing adapter removal_
4. **bowtie2.sh**: _index and alignment to the_ C. elegans _genome_
5. **sam_to_bam(1).sh**: _conversion of bowtie2 .sam output to .bam format_
6. **miRDeep2.sh**: _pre-processing and identification novel miRNAs_
7. **featureCounts_annotation(1).R**: _annotation and quantification of aligned reads_
8. **edgeR_miRNA_analysis.R**: _normalisation and analysis of miRNAs, compared against N2

### Polysome-seq

HISAT gene duplicates: **gene_duplicates.txt** 

1. **fastQC(2).sh** - _initial assessment of read quality_
2. **cutadapt(2).sh** - _removal of sequencing adapters_
3. **fastQC(2).sh** - _assessment of read quality folowing adapter removal_
4. **hisat2.sh** - _alignment of trimmed reads to the_ C. elegans _genome_
5. **remove_duplicates.sh** - _removal of gene duplicates generated during HISAT2 alignment_
6. **sam_to_bam(2).sh** - _conversion of HISAT2 .sam output to .bam format_
7. **filter_pairs.sh** - _filtering of processed reads to include matching pairs only_
8. **sort_bam.sh** - _sorting of .bam files for annotation_
9. **featureCounts_annotation(2).R** - _annotation and quantification of aligned reads_
10. **edgeR_polysome_analysis.R** - _normalisation and analysis of mRNAs_



## Citation
