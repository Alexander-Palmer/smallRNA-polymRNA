#Normalise and analyse processed small RNA-seq data
###################################################

#Programs
install.packages("edgeR")
library(edgeR)

#Input
a <- read.csv('meta_smallRNA_counts.csv', header = T)
b <- read.csv('meta_smallRNA_miRNAs.csv', header = T)

a <- as.data.frame(a)
a[is.na(a)] <- 0
sequences <- b$GeneID
levels <- letters[1:14]
levels[c(1,2)] <- "A"
levels[c(3,4)] <- "B"
levels[c(5,6)] <- "C"
levels[c(7,8)] <- "D"
levels[c(9,10)] <- "E"
levels[c(11,12)] <- "F"
levels[c(13,14)] <- "G

##############
#Processing_1#
##############

a <- DGEList(counts=a, group=levels, genes=sequences)

keep <- rowSums(cpm(a)>1) >= 1
a <- a[keep, , keep.lib.sizes=FALSE]

a <- calcNormFactors(a)
a <- estimateDisp(a)

design <- model.matrix(~0+levels, data=a$samples)
colnames(design) <- levels(a$samples$levels)
fit <- glmQLFit(a, design)

write.table(fit$genes, "meta_smallRNA_threshold_miRNAs.csv", sep=",", col.names = NA, quote = FALSE)
write.table(fit$fitted.values, "meta_smallRNA_threshold_normalised.csv", sep=",", col.names = NA, quote = FALSE)

##############
#Processing_2#
##############

#A-284_intestine_ALG1_vs_N2
#Positive logFC reflects positive enrichment in A-284_intestine_ALG1
Plcehldr <- glmQLFTest(fit, contrast=c(1,0,0,0,0,0,-1))
A284_N2_1 <- topTags(Plcehldr, n=Inf)
A284_N2_1 <- as.data.frame(A284_N2_1)
A284_N2_2 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
A284_N2_2 <- as.data.frame(A284_N2_2)

#A-285_intestine_ALG2_vs_N2
#Positive logFC reflects positive enrichment in A-285_intestine_ALG2
Plcehldr <- glmQLFTest(fit, contrast=c(0,1,0,0,0,0,-1))
A285_N2_1 <- topTags(Plcehldr, n=Inf)
A285_N2_1 <- as.data.frame(A285_N2_1)
A285_N2_2 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
A285_N2_2 <- as.data.frame(A285_N2_2)


#A-286_neuron_ALG1_vs_N2
#Positive logFC reflects positive enrichment in A-286_neuron_ALG1
Plcehldr <- glmQLFTest(fit, contrast=c(0,0,1,0,0,0,-1))
A286_N2_1 <- topTags(Plcehldr, n=Inf)
A286_N2_1 <- as.data.frame(A286_N2_1)
A286_N2_2 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
A286_N2_2 <- as.data.frame(A286_N2_2)

#A-287_neuron_ALG2_vs_N2
#Positive logFC reflects positive enrichment in A-287_neuron_ALG2
Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,1,0,0,-1))
A287_N2_1 <- topTags(Plcehldr, n=Inf)
A287_N2_1 <- as.data.frame(A287_N2_1)
A287_N2_2 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
A287_N2_2 <- as.data.frame(A287_N2_2)


#A-290_bwm_ALG1_vs_N2
#Positive logFC reflects positive enrichment in A-290_bwm_ALG1
Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,1,0,-1))
A290_N2_1 <- topTags(Plcehldr, n=Inf)
A290_N2_1 <- as.data.frame(A290_N2_1)
A290_N2_2 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
A290_N2_2 <- as.data.frame(A290_N2_2)

#A-291_bwm_ALG2_vs_N2
#Positive logFC reflects positive enrichment in A-291_bwm_ALG2
Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,0,1,-1))
A291_N2_1 <- topTags(Plcehldr, n=Inf)
A291_N2_1 <- as.data.frame(A291_N2_1)
A291_N2_2 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
A291_N2_2 <- as.data.frame(A291_N2_2)

##############
#Processing_3#
##############

sig.genes <- paste0("^", d$sig.genes, "$")

c <- read.csv("meta_smallRNA_threshold_normalised", header = TRUE)

sig.polysome.genes <- c[grep(paste(sig.genes, collapse = "|"), c[,1]), ]

##################
#Printing_results#
##################

write.table(A284_N2_2, "Intestine_ALG-1_vs_N2.csv", sep=",", col.names = NA, quote = FALSE)
write.table(A284_N2_1, "Intestine_ALG-1_vs_N2_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(A285_N2_2, "Intestine_ALG-2_vs_N2.csv", sep=",", col.names = NA, quote = FALSE)
write.table(A285_N2_1, "Intestine_ALG-2_vs_N2_all.csv", sep=",", col.names = NA, quote = FALSE)

write.table(A286_N2_2, "Neuron_ALG-1_vs_N2.csv", sep=",", col.names = NA, quote = FALSE)
write.table(A286_N2_1, "Neuron_ALG-1_vs_N2_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(A287_N2_2, "Neuron_ALG-2_vs_N2.csv", sep=",", col.names = NA, quote = FALSE)
write.table(A287_N2_1, "Neuron_ALG-2_vs_N2_all.csv", sep=",", col.names = NA, quote = FALSE)

write.table(A290_N2_2, "BWM_ALG-1_vs_N2.csv", sep=",", col.names = NA, quote = FALSE)
write.table(A290_N2_1, "BWM_ALG-1_vs_N2_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(A291_N2_2, "BWM_ALG-2_vs_N2.csv", sep=",", col.names = NA, quote = FALSE)
write.table(A291_N2_1, "BWM_ALG-2_vs_N2_all.csv", sep=",", col.names = NA, quote = FALSE)

write.table(c, "Sig_smallRNA_genes", sep=",", col.names = NA, quote = FALSE)
