#Annotates processed reads by C. elegans genome features
###########################################################################################################

#Download and apply requried packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")
library(Rsubread)

#Data input
proc.files <- c("S996-1-filtered-sorted.bam", "S996-2-filtered-sorted.bam", "S996-3-filtered-sorted.bam",
           "S997-1-filtered-sorted.bam", "S997-2-filtered-sorted.bam", "S997-3-filtered-sorted.bam",
           "S998-1-filtered-sorted.bam", "S998-2-filtered-sorted.bam", "S998-3-filtered-sorted.bam",
           "S999-1-filtered-sorted.bam", "S999-2-filtered-sorted.bam", "S999-3-filtered-sorted.bam",
           "S1000-1-filtered-sorted.bam", "S1000-2-filtered-sorted.bam", "S1000-3-filtered-sorted.bam",
           "S1045-1-filtered-sorted.bam", "S1045-2-filtered-sorted.bam", "S1045-3-filtered-sorted.bam",
           "S1046-1-filtered-sorted.bam", "S1046-2-filtered-sorted.bam", "S1046-3-filtered-sorted.bam",
           "S1047-1-filtered-sorted.bam", "S1047-2-filtered-sorted.bam", "S1047-3-filtered-sorted.bam")

#featureCounts processing
meta_results <-featureCounts(files=proc.files,
                       annot.ext="Caenorhabditis_elegans.WBcel235.96.gtf",
                       isGTFAnnotationFile=TRUE, GTF.featureType="gene",GTF.attrType="gene_id",
                       isPairedEnd=T, nthreads=12, useMetaFeatures=T,
                       genome="Caenorhabditis_elegans.WBcel235.dna.toplevel.fa",
                       juncCounts = T, allowMultiOverlap=T)

#Data output
write.table(x=data.frame(meta_results$annotation[,c("GeneID","Length")],
                         meta_results$counts,stringsAsFactors=FALSE),
                        file="counts_paired_final.txt",quote=FALSE,sep="\t")
