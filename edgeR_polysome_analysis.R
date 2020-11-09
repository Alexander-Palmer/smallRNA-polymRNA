###Quantify polysome-seq###

#Programs
install.packages("edgeR")
library(edgeR)

#Input
a <- read.csv('meta_polysome_hisat_clean.csv', header = T)
b <- read.csv('meta_polysome_hisat_identity.csv', header = T)
c <- read.csv('meta_polysome_normalised_no_thresh.csv', header = T)

c[paste("Av_996")] <- ((c[["S996.1"]] + c[["S996.2"]] + c[["S996.3"]]) / 3)
c[paste("Av_997")] <- ((c[["S997.1"]] + c[["S997.2"]] + c[["S997.3"]]) / 3)
c[paste("Av_998")] <- ((c[["S998.1"]] + c[["S998.2"]] + c[["S998.3"]]) / 3)
c[paste("Av_999")] <- ((c[["S999.1"]] + c[["S999.2"]] + c[["S999.3"]]) / 3)
c[paste("Av_1000")] <- ((c[["S1000.1"]] + c[["S1000.2"]] + c[["S1000.3"]]) / 3)
c[paste("Av_1045")] <- ((c[["S1045.1"]] + c[["S1045.2"]] + c[["S1045.3"]]) / 3)
c[paste("Av_1046")] <- ((c[["S1046.1"]] + c[["S1045.2"]] + c[["S1046.3"]]) / 3)
c[paste("Av_1047")] <- ((c[["S1047.1"]] + c[["S1045.2"]] + c[["S1047.3"]]) / 3)
c[paste("Av_combined")] <- ((c[["Av_996"]] + c[["Av_997"]] + c[["Av_998"]] + c[["Av_999"]] + c[["Av_1000"]]
                             + c[["Av_1045"]] + c[["Av_1046"]] + c[["Av_1047"]]) / 8)

write.table(c, "meta_polysome_normalised_no_thresh_averages.csv", sep=",", col.names = NA, quote = FALSE)

a <- as.data.frame(a)
a[is.na(a)] <- 0
sequences <- b$GeneID
levels <- letters[1:24]
levels[c(1,2,3)] <- "A"
levels[c(4,5,6)] <- "B"
levels[c(7,8,9)] <- "C"
levels[c(10,11,12)] <- "D"
levels[c(13,14,15)] <- "E"
levels[c(16,17,18)] <- "F"
levels[c(19,20,21)] <- "G"
levels[c(22,23,24)] <- "H"

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

write.table(fit$genes, "meta_polysome_normalised_names.csv", sep=",", col.names = NA, quote = FALSE)
write.table(fit$fitted.values, "meta_polysome_normalised.csv", sep=",", col.names = NA, quote = FALSE)

##############
#Processing_2#
##############

#S996
{
  #S996_vs_S997
  Plcehldr <- glmQLFTest(fit, contrast=c(1,-1,0,0,0,0,0,0))
  a1 <- topTags(Plcehldr, n=Inf)
  a1 <- as.data.frame(a1)
  b1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  b1 <- as.data.frame(b1)

  #S996_vs_S998
  Plcehldr <- glmQLFTest(fit, contrast=c(1,0,-1,0,0,0,0,0))
  c1 <- topTags(Plcehldr, n=Inf)
  c1 <- as.data.frame(c1)
  d1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  d1 <- as.data.frame(d1)

  #S996_vs_S999
  Plcehldr <- glmQLFTest(fit, contrast=c(1,0,0,-1,0,0,0,0))
  e1 <- topTags(Plcehldr, n=Inf)
  e1 <- as.data.frame(e1)
  f1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  f1 <- as.data.frame(f1)

  #S996_vs_S1000
  Plcehldr <- glmQLFTest(fit, contrast=c(1,0,0,0,-1,0,0,0))
  g1 <- topTags(Plcehldr, n=Inf)
  g1 <- as.data.frame(g1)
  h1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  h1 <- as.data.frame(h1)

  #S996_vs_S1045
  Plcehldr <- glmQLFTest(fit, contrast=c(1,0,0,0,0,-1,0,0))
  i1 <- topTags(Plcehldr, n=Inf)
  i1 <- as.data.frame(i1)
  j1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  j1 <- as.data.frame(j1)

  #S996_vs_S1046
  Plcehldr <- glmQLFTest(fit, contrast=c(1,0,0,0,0,0,-1,0))
  k1 <- topTags(Plcehldr, n=Inf)
  k1 <- as.data.frame(k1)
  l1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  l1 <- as.data.frame(l1)

  #S996_vs_S1047
  Plcehldr <- glmQLFTest(fit, contrast=c(1,0,0,0,0,0,0,-1))
  m1 <- topTags(Plcehldr, n=Inf)
  m1 <- as.data.frame(m1)
  n1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  n1 <- as.data.frame(n1)
}


#S997
{
  #S997_vs_S996
  Plcehldr <- glmQLFTest(fit, contrast=c(-1,1,0,0,0,0,0,0))
  o1 <- topTags(Plcehldr, n=Inf)
  o1 <- as.data.frame(o1)
  p1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  p1 <- as.data.frame(p1)

  #S997_vs_S998
  Plcehldr <- glmQLFTest(fit, contrast=c(0,1,-1,0,0,0,0,0))
  q1 <- topTags(Plcehldr, n=Inf)
  q1 <- as.data.frame(q1)
  r1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  r1 <- as.data.frame(r1)

  #S997_vs_S999
  Plcehldr <- glmQLFTest(fit, contrast=c(0,1,0,-1,0,0,0,0))
  s1 <- topTags(Plcehldr, n=Inf)
  s1 <- as.data.frame(s1)
  t1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  t1 <- as.data.frame(t1)

  #S997_vs_S1000
  Plcehldr <- glmQLFTest(fit, contrast=c(0,1,0,0,-1,0,0,0))
  u1 <- topTags(Plcehldr, n=Inf)
  u1 <- as.data.frame(u1)
  v1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  v1 <- as.data.frame(v1)

  #S997_vs_S1045
  Plcehldr <- glmQLFTest(fit, contrast=c(0,1,0,0,0,-1,0,0))
  w1 <- topTags(Plcehldr, n=Inf)
  w1 <- as.data.frame(w1)
  x1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  x1 <- as.data.frame(x1)

  #S997_vs_S1046
  Plcehldr <- glmQLFTest(fit, contrast=c(0,1,0,0,0,0,-1,0))
  y1 <- topTags(Plcehldr, n=Inf)
  y1 <- as.data.frame(y1)
  z1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  z1 <- as.data.frame(z1)

  #S997_vs_S1047
  Plcehldr <- glmQLFTest(fit, contrast=c(0,1,0,0,0,0,0,-1))
  aa1 <- topTags(Plcehldr, n=Inf)
  aa1 <- as.data.frame(aa1)
  ab1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  ab1 <- as.data.frame(ab1)
}


#S998
{
  #S998_vs_S996
  Plcehldr <- glmQLFTest(fit, contrast=c(-1,0,1,0,0,0,0,0))
  ac1 <- topTags(Plcehldr, n=Inf)
  ac1 <- as.data.frame(ac1)
  ad1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  ad1 <- as.data.frame(ad1)

  #S998_vs_S997
  Plcehldr <- glmQLFTest(fit, contrast=c(0,-1,1,0,0,0,0,0))
  ae1 <- topTags(Plcehldr, n=Inf)
  ae1 <- as.data.frame(ae1)
  af1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  af1 <- as.data.frame(af1)

  #S998_vs_S999
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,1,-1,0,0,0,0))
  ag1 <- topTags(Plcehldr, n=Inf)
  ag1 <- as.data.frame(ag1)
  ah1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  ah1 <- as.data.frame(ah1)

  #S998_vs_S1000
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,1,0,-1,0,0,0))
  ai1 <- topTags(Plcehldr, n=Inf)
  ai1 <- as.data.frame(ai1)
  aj1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  aj1 <- as.data.frame(aj1)

  #S998_vs_S1045
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,1,0,0,-1,0,0))
  ak1 <- topTags(Plcehldr, n=Inf)
  ak1 <- as.data.frame(ak1)
  al1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  al1 <- as.data.frame(al1)

  #S998_vs_S1046
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,1,0,0,0,-1,0))
  am1 <- topTags(Plcehldr, n=Inf)
  am1 <- as.data.frame(am1)
  an1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  an1 <- as.data.frame(an1)

  #S998_vs_S1047
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,1,0,0,0,0,-1))
  ao1 <- topTags(Plcehldr, n=Inf)
  ao1 <- as.data.frame(ao1)
  ap1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  ap1 <- as.data.frame(ap1)
}


#S999
{
  #S999_vs_S996
  Plcehldr <- glmQLFTest(fit, contrast=c(-1,0,0,1,0,0,0,0))
  aq1 <- topTags(Plcehldr, n=Inf)
  aq1 <- as.data.frame(aq1)
  ar1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  ar1 <- as.data.frame(ar1)

  #S999_vs_S997
  Plcehldr <- glmQLFTest(fit, contrast=c(0,-1,0,1,0,0,0,0))
  as1 <- topTags(Plcehldr, n=Inf)
  as1 <- as.data.frame(as1)
  at1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  at1 <- as.data.frame(at1)

  #S999_vs_S998
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,-1,1,0,0,0,0))
  au1 <- topTags(Plcehldr, n=Inf)
  au1 <- as.data.frame(au1)
  av1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  av1 <- as.data.frame(av1)

  #S999_vs_S1000
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,1,-1,0,0,0))
  aw1 <- topTags(Plcehldr, n=Inf)
  aw1 <- as.data.frame(aw1)
  ax1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  ax1 <- as.data.frame(ax1)

  #S999_vs_S1045
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,1,0,-1,0,0))
  ay1 <- topTags(Plcehldr, n=Inf)
  ay2 <- as.data.frame(ay1)
  az1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  az1 <- as.data.frame(az1)

  #S999_vs_S1046
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,1,0,0,-1,0))
  ba1 <- topTags(Plcehldr, n=Inf)
  ba2 <- as.data.frame(ba1)
  bb1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  bb1 <- as.data.frame(bb1)

  #S999_vs_S1047
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,1,0,0,0,-1))
  bc1 <- topTags(Plcehldr, n=Inf)
  bc2 <- as.data.frame(bc1)
  bd1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  bd1 <- as.data.frame(bd1)
}


#S1000
{
  #S1000_vs_S996
  Plcehldr <- glmQLFTest(fit, contrast=c(-1,0,0,0,1,0,0,0))
  be1 <- topTags(Plcehldr, n=Inf)
  be1 <- as.data.frame(be1)
  bf1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  bf1 <- as.data.frame(bf1)

  #S1000_vs_S997
  Plcehldr <- glmQLFTest(fit, contrast=c(0,-1,0,0,1,0,0,0))
  bg1 <- topTags(Plcehldr, n=Inf)
  bg1 <- as.data.frame(bg1)
  bh1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  bh1 <- as.data.frame(bh1)

  #S1000_vs_S998
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,-1,0,1,0,0,0))
  bi1 <- topTags(Plcehldr, n=Inf)
  bi1 <- as.data.frame(bi1)
  bj1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  bj1 <- as.data.frame(bj1)

  #S1000_vs_S999
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,-1,1,0,0,0))
  bk1 <- topTags(Plcehldr, n=Inf)
  bk1 <- as.data.frame(bk1)
  bl1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  bl1 <- as.data.frame(bl1)

  #S1000_vs_S1045
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,1,-1,0,0))
  bm1 <- topTags(Plcehldr, n=Inf)
  bm1 <- as.data.frame(bm1)
  bn1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  bn1 <- as.data.frame(bn1)

  #S1000_vs_S1046
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,1,0,-1,0))
  bo1 <- topTags(Plcehldr, n=Inf)
  bo1 <- as.data.frame(bo1)
  bp1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  bp1 <- as.data.frame(bp1)

  #S1000_vs_S1047
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,1,0,0,-1))
  bq1 <- topTags(Plcehldr, n=Inf)
  bq1 <- as.data.frame(bq1)
  br1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  br1 <- as.data.frame(br1)
}


#S1045
{
  #S1045_vs_S996
  Plcehldr <- glmQLFTest(fit, contrast=c(-1,0,0,0,0,1,0,0))
  bs1 <- topTags(Plcehldr, n=Inf)
  bs1 <- as.data.frame(bs1)
  bt1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  bt1 <- as.data.frame(bt1)

  #S1045_vs_S997
  Plcehldr <- glmQLFTest(fit, contrast=c(0,-1,0,0,0,1,0,0))
  bu1 <- topTags(Plcehldr, n=Inf)
  bu1 <- as.data.frame(bu1)
  bv1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  bv1 <- as.data.frame(bv1)

  #S1045_vs_S998
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,-1,0,0,1,0,0))
  bw1 <- topTags(Plcehldr, n=Inf)
  bw1 <- as.data.frame(bw1)
  bx1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  bx1 <- as.data.frame(bx1)

  #S1045_vs_S999
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,-1,0,1,0,0))
  by1 <- topTags(Plcehldr, n=Inf)
  by1 <- as.data.frame(by1)
  bz1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  bz1 <- as.data.frame(bz1)

  #S1045_vs_S1000
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,-1,1,0,0))
  ca1 <- topTags(Plcehldr, n=Inf)
  ca1 <- as.data.frame(ca1)
  cb1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  cb1 <- as.data.frame(cb1)

  #S1045_vs_S1046
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,0,1,-1,0))
  cc1 <- topTags(Plcehldr, n=Inf)
  cc1 <- as.data.frame(cc1)
  cd1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  cd1 <- as.data.frame(cd1)

  #S1045_vs_S1047
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,0,1,0,-1))
  ce1 <- topTags(Plcehldr, n=Inf)
  ce1 <- as.data.frame(ce1)
  cf1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  cf1 <- as.data.frame(cf1)
}


#S1046
{
  #S1045_vs_S996
  Plcehldr <- glmQLFTest(fit, contrast=c(-1,0,0,0,0,0,1,0))
  cg1 <- topTags(Plcehldr, n=Inf)
  cg1 <- as.data.frame(cg1)
  ch1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  ch1 <- as.data.frame(ch1)

  #S1046_vs_S997
  Plcehldr <- glmQLFTest(fit, contrast=c(0,-1,0,0,0,0,1,0))
  ci1 <- topTags(Plcehldr, n=Inf)
  ci1 <- as.data.frame(ci1)
  cj1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  cj1 <- as.data.frame(cj1)

  #S1046_vs_S998
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,-1,0,0,0,1,0))
  ck1 <- topTags(Plcehldr, n=Inf)
  ck1 <- as.data.frame(ck1)
  cl1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  cl1 <- as.data.frame(cl1)

  #S1046_vs_S999
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,-1,0,0,1,0))
  cm1 <- topTags(Plcehldr, n=Inf)
  cm1 <- as.data.frame(cm1)
  cn1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  cn1 <- as.data.frame(cn1)

  #S1046_vs_S1000
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,-1,0,1,0))
  co1 <- topTags(Plcehldr, n=Inf)
  co1 <- as.data.frame(co1)
  cp1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  cp1 <- as.data.frame(cp1)

  #S1046_vs_S1045
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,0,-1,1,0))
  cq1 <- topTags(Plcehldr, n=Inf)
  cq1 <- as.data.frame(cq1)
  cr1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  cr1 <- as.data.frame(cr1)

  #S1046_vs_S1047
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,0,0,1,-1))
  cs1 <- topTags(Plcehldr, n=Inf)
  cs1 <- as.data.frame(cs1)
  ct1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  ct1 <- as.data.frame(ct1)
}


#S1047
{
  #S1047_vs_S996
  Plcehldr <- glmQLFTest(fit, contrast=c(-1,0,0,0,0,0,0,1))
  cu1 <- topTags(Plcehldr, n=Inf)
  cu1 <- as.data.frame(cu1)
  cv1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  cv1 <- as.data.frame(cv1)

  #S1047_vs_S997
  Plcehldr <- glmQLFTest(fit, contrast=c(0,-1,0,0,0,0,0,1))
  cw1 <- topTags(Plcehldr, n=Inf)
  cw1 <- as.data.frame(cw1)
  cx1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  cx1 <- as.data.frame(cx1)

  #S1047_vs_S998
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,-1,0,0,0,0,1))
  cy1 <- topTags(Plcehldr, n=Inf)
  cy1 <- as.data.frame(cy1)
  cz1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  cz1 <- as.data.frame(cz1)

  #S1047_vs_S999
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,-1,0,0,0,1))
  da1 <- topTags(Plcehldr, n=Inf)
  da1 <- as.data.frame(da1)
  db1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  db1 <- as.data.frame(db1)

  #S1047_vs_S1000
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,-1,0,0,1))
  dc1 <- topTags(Plcehldr, n=Inf)
  dc1 <- as.data.frame(dc1)
  dd1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  dd1 <- as.data.frame(dd1)

  #S1047_vs_S1045
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,0,-1,0,1))
  de1 <- topTags(Plcehldr, n=Inf)
  de1 <- as.data.frame(de1)
  df1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  df1 <- as.data.frame(df1)

  #S1047_vs_S1046
  Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,0,0,-1,1))
  dg1 <- topTags(Plcehldr, n=Inf)
  dg1 <- as.data.frame(dg1)
  dh1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
  dh1 <- as.data.frame(dh1)
}

##############
#Processing_3#
##############

sig.genes <- Reduce(union, list(b1$genes, d1$genes, f1$genes, h1$genes, j1$genes, l1$genes, n1$genes, p1$genes, r1$genes,
                                t1$genes, v1$genes, x1$genes, z1$genes, ab1$genes, ad1$genes, af1$genes, ah1$genes,
                                aj1$genes, al1$genes, an1$genes, ap1$genes, ar1$genes, at1$genes, av1$genes, ax1$genes,
                                az1$genes, bb1$genes, bd1$genes, bf1$genes, bh1$genes, bj1$genes, bl1$genes, bn1$genes,
                                bp1$genes, br1$genes, bt1$genes, bv1$genes, bx1$genes, bz1$genes, cb1$genes, cd1$genes,
                                cf1$genes, ch1$genes, cj1$genes, cl1$genes, cn1$genes, cp1$genes, cr1$genes, ct1$genes,
                                cv1$genes, cx1$genes, cz1$genes, db1$genes, dd1$genes, df1$genes, dh1$genes))
sig.genes <- paste0("^", sig.genes, "$")

d <- read.csv("meta_polysome_normalised_no_thresh_averages.csv", header = TRUE)

sig.polysome.genes <- d[grep(paste(sig.genes, collapse = "|"), d[,1]), ]

##################
#Printing_results#
##################

write.table(b1, "996x997.csv", sep=",", col.names = NA, quote = FALSE)
write.table(d1, "996x998.csv", sep=",", col.names = NA, quote = FALSE)
write.table(f1, "996x999.csv", sep=",", col.names = NA, quote = FALSE)
write.table(h1, "996x1000.csv", sep=",", col.names = NA, quote = FALSE)
write.table(j1, "996x1045.csv", sep=",", col.names = NA, quote = FALSE)
write.table(l1, "996x1046.csv", sep=",", col.names = NA, quote = FALSE)
write.table(n1, "996x1047.csv", sep=",", col.names = NA, quote = FALSE)
write.table(p1, "997x996.csv", sep=",", col.names = NA, quote = FALSE)
write.table(r1, "997x998.csv", sep=",", col.names = NA, quote = FALSE)
write.table(t1, "997x999.csv", sep=",", col.names = NA, quote = FALSE)
write.table(v1, "997x1000.csv", sep=",", col.names = NA, quote = FALSE)
write.table(x1, "997x1045.csv", sep=",", col.names = NA, quote = FALSE)
write.table(z1, "997x1046.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ab1, "997x1047.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ad1, "998x996.csv", sep=",", col.names = NA, quote = FALSE)
write.table(af1, "998x997.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ah1, "998x999.csv", sep=",", col.names = NA, quote = FALSE)
write.table(aj1, "998x1000.csv", sep=",", col.names = NA, quote = FALSE)
write.table(al1, "998x1045.csv", sep=",", col.names = NA, quote = FALSE)
write.table(an1, "998x1046.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ap1, "998x1047.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ar1, "999x996.csv", sep=",", col.names = NA, quote = FALSE)
write.table(at1, "999x997.csv", sep=",", col.names = NA, quote = FALSE)
write.table(av1, "999x998.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ax1, "999x1000.csv", sep=",", col.names = NA, quote = FALSE)
write.table(az1, "999x1045.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bb1, "999x1046.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bd1, "999x1047.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bf1, "1000x996.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bh1, "1000x997.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bj1, "1000x998.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bl1, "1000x999.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bn1, "1000x1045.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bp1, "1000x1046.csv", sep=",", col.names = NA, quote = FALSE)
write.table(br1, "1000x1047.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bt1, "1045x996.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bv1, "1045x997.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bx1, "1045x998.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bz1, "1045x999.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cb1, "1045x1000.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cd1, "1045x1046.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cf1, "1045x1047.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ch1, "1046x996.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cj1, "1046x997.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cl1, "1046x998.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cn1, "1046x999.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cp1, "1046x1000.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cr1, "1046x1045.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ct1, "1046x1047.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cv1, "1047x996.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cx1, "1047x997.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cz1, "1047x998.csv", sep=",", col.names = NA, quote = FALSE)
write.table(db1, "1047x999.csv", sep=",", col.names = NA, quote = FALSE)
write.table(dd1, "1047x1000.csv", sep=",", col.names = NA, quote = FALSE)
write.table(df1, "1047x1045.csv", sep=",", col.names = NA, quote = FALSE)
write.table(dh1, "1047x1046.csv", sep=",", col.names = NA, quote = FALSE)

###

write.table(a1, "996x997_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(c1, "996x998_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(e1, "996x999_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(g1, "996x1000_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(i1, "996x1045_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(k1, "996x1046_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(m1, "996x1047_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(o1, "997x996_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(q1, "997x998_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(s1, "997x999_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(u1, "997x1000_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(w1, "997x1045_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(y1, "997x1046_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(aa1, "997x1047_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ac1, "998x996_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ae1, "998x997_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ag1, "998x999_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ai1, "998x1000_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ak1, "998x1045_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(am1, "998x1046_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ao1, "998x1047_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(aq1, "999x996_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(as1, "999x997_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(au1, "999x998_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(aw1, "999x1000_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ay1, "999x1045_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ba1, "999x1046_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bc1, "999x1047_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(be1, "1000x996_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bg1, "1000x997_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bi1, "1000x998_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bk1, "1000x999_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bm1, "1000x1045_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bo1, "1000x1046_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bq1, "1000x1047_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bs1, "1045x996_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bu1, "1045x997_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(bw1, "1045x998_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(by1, "1045x999_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ca1, "1045x1000_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cc1, "1045x1046_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ce1, "1045x1047_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cg1, "1046x996_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ci1, "1046x997_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(ck1, "1046x998_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cm1, "1046x999_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(co1, "1046x1000_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cq1, "1046x1045_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cs1, "1046x1047_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cu1, "1047x996_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cw1, "1047x997_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(cy1, "1047x998_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(da1, "1047x999_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(dc1, "1047x1000_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(de1, "1047x1045_all.csv", sep=",", col.names = NA, quote = FALSE)
write.table(dg1, "1047x1046_all.csv", sep=",", col.names = NA, quote = FALSE)








#sig.polysome.genes <- sig.polysome.genes[!grepl("C03B1.10", sig.polysome.genes$target_id),]



#Remove individual columns
{
sig.polysome.genes$S996.1.tpm <- NULL
sig.polysome.genes$S996.2.tpm <- NULL
sig.polysome.genes$S996.3.tpm <- NULL

sig.polysome.genes$S997.1.tpm <- NULL
sig.polysome.genes$S997.2.tpm <- NULL
sig.polysome.genes$S997.3.tpm <- NULL

sig.polysome.genes$S998.1.tpm <- NULL
sig.polysome.genes$S998.2.tpm <- NULL
sig.polysome.genes$S998.3.tpm <- NULL

sig.polysome.genes$S999.1.tpm <- NULL
sig.polysome.genes$S999.2.tpm <- NULL
sig.polysome.genes$S999.3.tpm <- NULL

sig.polysome.genes$S1000.1.tpm <- NULL
sig.polysome.genes$S1000.2.tpm <- NULL
sig.polysome.genes$S1000.3.tpm <- NULL

sig.polysome.genes$S1045.1.tpm <- NULL
sig.polysome.genes$S1045.2.tpm <- NULL
sig.polysome.genes$S1045.3.tpm <- NULL

sig.polysome.genes$S1046.1.tpm <- NULL
sig.polysome.genes$S1046.2.tpm <- NULL
sig.polysome.genes$S1046.3.tpm <- NULL

sig.polysome.genes$S1047.1.tpm <- NULL
sig.polysome.genes$S1047.2.tpm <- NULL
sig.polysome.genes$S1047.3.tpm <- NULL
}

sig.polysome.genes$Av_996 <- log2(sig.polysome.genes$Av_996)
sig.polysome.genes$Av_997 <- log2(sig.polysome.genes$Av_997)
sig.polysome.genes$Av_999 <- log2(sig.polysome.genes$Av_999)
sig.polysome.genes$Av_1045 <- log2(sig.polysome.genes$Av_1045)
sig.polysome.genes$Av_1046 <- log2(sig.polysome.genes$Av_1046)
sig.polysome.genes$Av_1047 <- log2(sig.polysome.genes$Av_1047)

#OR
{
sig.polysome.genes$S996.1 <- log2(sig.polysome.genes$S996.1)
sig.polysome.genes$S996.2 <- log2(sig.polysome.genes$S996.2)
sig.polysome.genes$S996.3 <- log2(sig.polysome.genes$S996.3)

sig.polysome.genes$S997.1 <- log2(sig.polysome.genes$S997.1)
sig.polysome.genes$S997.2 <- log2(sig.polysome.genes$S997.2)
sig.polysome.genes$S997.3 <- log2(sig.polysome.genes$S997.3)

sig.polysome.genes$S999.1 <- log2(sig.polysome.genes$S999.1)
sig.polysome.genes$S999.2 <- log2(sig.polysome.genes$S999.2)
sig.polysome.genes$S999.3 <- log2(sig.polysome.genes$S999.3)

sig.polysome.genes$S1045.1 <- log2(sig.polysome.genes$S1045.1)
sig.polysome.genes$S1045.2 <- log2(sig.polysome.genes$S1045.2)
sig.polysome.genes$S1045.3 <- log2(sig.polysome.genes$S1045.3)

sig.polysome.genes$S1046.1 <- log2(sig.polysome.genes$S1046.1)
sig.polysome.genes$S1046.2 <- log2(sig.polysome.genes$S1046.2)
sig.polysome.genes$S1046.3 <- log2(sig.polysome.genes$S1046.3)

sig.polysome.genes$S1047.1 <- log2(sig.polysome.genes$S1047.1)
sig.polysome.genes$S1047.2 <- log2(sig.polysome.genes$S1047.2)
sig.polysome.genes$S1047.3 <- log2(sig.polysome.genes$S1047.3)
}

sig.polysome.genes[sig.polysome.genes < 0] <- NA
sig.polysome.genes[is.na(sig.polysome.genes)] <- 0

sig.polysome.genes$target_id_mRNA <- re.sub(r'^((?:[^.]*\.){1}[^.]*)\..*', r'\1', sig.polysome.genes$target_id)
