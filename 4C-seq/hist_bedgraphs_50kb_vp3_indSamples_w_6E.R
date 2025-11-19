# Usage:
#   Run in a directory containing gzipped 50 kb bedGraph files:
#     S9_6E_1_EcoRI_3.50kb.bedGraph.gz, ...
#
#   This script:
#     - pools all samples to estimate a global mean,
#     - computes Poisson p-values per bin (ppois, lower.tail=FALSE),
#     - defines significant bins as -log10(p) > 5,
#     - merges adjacent significant bins into regions >= 50 kb,
#     - writes:
#         <sample>.tem.50kb.bed
#         <sample>.tem.50kb.ppois.bedgraph

library(GenomicRanges)
library(ggplot2)

a <- rtracklayer::import("../BedGraph/S12_9E_1_EcoRI_3.50kb.bedGraph.gz")
b <- rtracklayer::import("../BedGraph/S13_9E_2_EcoRI_3.50kb.bedGraph.gz")
d <- rtracklayer::import("../BedGraph/S14_9E_3_EcoRI_3.50kb.bedGraph.gz")
e <- rtracklayer::import("../BedGraph/S15_9E_4_EcoRI_3.50kb.bedGraph.gz")
f <- rtracklayer::import("../BedGraph/S9_6E_1_EcoRI_3.50kb.bedGraph.gz")
g <- rtracklayer::import("../BedGraph/S10_6E_2_EcoRI_3.50kb.bedGraph.gz")
h <- rtracklayer::import("../BedGraph/S11_6E_3_EcoRI_3.50kb.bedGraph.gz")

f <- dropSeqlevels(f, "chr4",pruning.mode="coarse")
g <- dropSeqlevels(g, "chr4",pruning.mode="coarse")
h <- dropSeqlevels(h, "chr4",pruning.mode="coarse")

######
alldat <- data.frame(Sample=c(rep("S12",length(a)),
                              rep("S13",length(b)),
                              rep("S14",length(d)),
                              rep("S15",length(e)),
                              rep("S9",length(f)),
                              rep("S10",length(g)),
                              rep("S11",length(h))),
                     meanNormReads=c(a$score,b$score,d$score,e$score,
                                     f$score,g$score,h$score))

meandat <- mean(alldat$meanNormReads) # 22.662

alldat$ppois <- ppois(alldat$meanNormReads,lambda=meandat,lower.tail=F)
alldat$score <- -log10(alldat$ppois)
alldat$score[which(alldat$score == Inf)] <- 350
alldat$sig <- ifelse(alldat$score > 5, TRUE, FALSE)

table(alldat$Sample,alldat$sig)
#FALSE  TRUE
#S10 19417  3181
#S11 25633  3194
#S12 33954  7105
#S13 24603  5960
#S14 50538  4443
#S15 26323  6609
#S9  20500  3104

ggplot(alldat, aes(meanNormReads, after_stat(density), fill = Sample)) +
     geom_histogram(binwidth = 10) + theme_bw() +facet_wrap(~Sample)

pdf("BySample_50kb_vp3_histograms.pdf")
ggplot(alldat, aes(meanNormReads, fill = sig)) +
  geom_histogram(binwidth = 10) + theme_bw() +facet_wrap(~Sample)
dev.off()

a$score <- alldat$score[which(alldat$Sample == "S12")]
a2 <- a[which(a$score > 5)]
mcols(a2) <- NULL
a2 <- reduce(a2, min.gapwidth=50001)
a2 <- a2[which(width(a2) > 50000)]
rtracklayer::export(a2,"S12.tem.50kb.bed")
rtracklayer::export(a,"S12.tem.50kb.ppois.bedgraph")
# 1370

b$score <- alldat$score[which(alldat$Sample == "S13")]
b2 <- b[which(b$score > 5)]
mcols(b2) <- NULL
b2 <- reduce(b2, min.gapwidth=50001)
b2 <- b2[which(width(b2) > 50000)]
rtracklayer::export(b2,"S13.tem.50kb.bed")
rtracklayer::export(b,"S13.tem.50kb.ppois.bedgraph")
# 1061

d$score <- alldat$score[which(alldat$Sample == "S14")]
d2 <- d[which(d$score > 5)]
mcols(d2) <- NULL
d2 <- reduce(d2, min.gapwidth=50001)
d2 <- d2[which(width(d2) > 50000)]
rtracklayer::export(d2,"S14.tem.50kb.bed")
rtracklayer::export(d,"S14.tem.50kb.ppois.bedgraph")
# 876

e$score <- alldat$score[which(alldat$Sample == "S15")]
e2 <- e[which(e$score > 5)]
mcols(e2) <- NULL
e2 <- reduce(e2, min.gapwidth=50001)
e2 <- e2[which(width(e2) > 50000)]
rtracklayer::export(e2,"S15.tem.50kb.bed")
rtracklayer::export(e,"S15.tem.50kb.ppois.bedgraph")
# 1204

f$score <- alldat$score[which(alldat$Sample == "S9")]
f2 <- f[which(f$score > 5)]
mcols(f2) <- NULL
f2 <- reduce(f2, min.gapwidth=50001)
f2 <- f2[which(width(f2) > 50000)]
rtracklayer::export(f2,"S9.tem.50kb.bed")
rtracklayer::export(f,"S9.tem.50kb.ppois.bedgraph")
# 473

g$score <- alldat$score[which(alldat$Sample == "S10")]
g2 <- g[which(g$score > 5)]
mcols(g2) <- NULL
g2 <- reduce(g2, min.gapwidth=50001)
g2 <- g2[which(width(g2) > 50000)]
rtracklayer::export(g2,"S10.tem.50kb.bed")
rtracklayer::export(g,"S10.tem.50kb.ppois.bedgraph")
# 465

h$score <- alldat$score[which(alldat$Sample == "S11")]
h2 <- h[which(h$score > 5)]
mcols(h2) <- NULL
h2 <- reduce(h2, min.gapwidth=50001)
h2 <- h2[which(width(h2) > 50000)]
rtracklayer::export(h2,"S11.tem.50kb.bed")
rtracklayer::export(h,"S11.tem.50kb.ppois.bedgraph")
# 470

