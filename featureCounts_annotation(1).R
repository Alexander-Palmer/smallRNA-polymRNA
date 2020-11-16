#Annotates small RNA-seq processed reads by C. elegans genome features
######################################################################

#Download and apply requried packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")
library(Rsubread)

#Data input
proc.files <- c("A-284-1_intestine_ALG1_clip.bam", "A-284-2_intestine_ALG1_clip.bam", 
                "A-285-1_intestine_ALG2_clip.bam", "A-285-2_intestine_ALG2_clip.bam",
                "A-286-1_neuron_ALG1_clip.bam", "A-286-2_neuron_ALG1_clip.bam",
                "A-287-1_neuron_ALG2_clip.bam", "A-287-2_neuron_ALG2_clip.bam",
                "A-290-1_bwm_ALG1_clip.bam", "A-290-2_bwm_ALG1_clip.bam",
                "A-291-1_bwm_ALG2_clip.bam", "A-291-2_bwm_ALG2_clip.bam")
                           
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
                        file="miRNA_counts_paired_final.txt",quote=FALSE,sep="\t")
