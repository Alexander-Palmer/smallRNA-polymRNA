#Normalise and analyse processed polysome-seq data
##################################################

#Programs
install.packages("edgeR")
library(edgeR)

#Input
a <- read.csv('meta_polysome_counts.csv', header = T)
b <- read.csv('meta_polysome_mRNAs.csv', header = T)

a <- as.data.frame(a)
a[is.na(a)] <- 0
sequences <- b$GeneID
levels <- letters[1:6]
levels[c(1,2,3)] <- "A"
levels[c(4,5,6)] <- "B"

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

write.table(fit$genes, "meta_polysome_threshold_mRNAs.csv", sep=",", col.names = NA, quote = FALSE)
write.table(fit$fitted.values, "meta_polysome_threshold_normalised.csv", sep=",", col.names = NA, quote = FALSE)

##############
#Processing_2#
##############

#S998_vs_S1000
#Positive logFC reflects upregulation in S998
Plcehldr <- glmQLFTest(fit, contrast=c(1,-1))
c <- topTags(Plcehldr, n=Inf)
c <- as.data.frame(c)
d <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
d <- as.data.frame(d)


##############
#Processing_3#
##############

sig.genes <- paste0("^", d$sig.genes, "$")

e <- read.csv("meta_polysome_threshold_normalised", header = TRUE)

sig.polysome.genes <- e[grep(paste(sig.genes, collapse = "|"), e[,1]), ]

##################
#Printing_results#
##################

write.table(d, "998x1000.csv", sep=",", col.names = NA, quote = FALSE)
write.table(c, "998x1000_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(e, "Sig_polysome_genes", sep=",", col.names = NA, quote = FALSE)
