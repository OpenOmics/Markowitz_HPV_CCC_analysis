# Usage:
#   Run in a directory containing per-sample 50 kb bedGraph files
#   named like:
#     <exp>_<cond>_<rep>_<digest>_<viewpoint>.50kb.bedGraph.gz
#
#   This script groups samples by experiment + viewpoint and writes:
#     <condition>_viewpoint<VP>.50kb.mean.bedgraph
#
# Requirements:
#   - R + packages: GenomicRanges, BSgenome.Hsapiens.UCSC.hg38, rtracklayer
#   - Input bedGraphs should have a 'score' column with normalized read counts.

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

inFiles <- list.files(pattern = "bedGraph.gz")

sampleNames <- gsub(".50kb.bedGraph.gz","",inFiles)
sampleInfo <- matrix(unlist(strsplit(sampleNames,"_")),ncol=5,byrow = T)
groupInfo <- paste0(sampleInfo[,2],"_viewpoint",sampleInfo[,5])
groups <- unique(groupInfo)

for (i in 1:length(groups)) {
  samps <- inFiles[which(groupInfo == groups[i])]
  sampData <- lapply(samps, rtracklayer::import)
  
  # This is to deal with the fact that we removed the empty rows in the original
  # bedgraphs and the objects no longer have equal numbers of rows
  tmp <- seqinfo(Hsapiens)
  seqlevels(tmp) <- seqlevels(sampData[[1]])[1:24]
  tiles <- tileGenome(tmp, tilewidth=5e4,
                      cut.last.tile.in.chrom=TRUE)
  tiles2 <- tiles
  mcols(tiles2) <- matrix(nrow=length(tiles2),ncol=length(samps),data=0)
  for (j in 1:length(samps)) {
        ov <- findOverlaps(tiles2, sampData[[j]])
        mcols(tiles2)[queryHits(ov),j] <- sampData[[j]]$score
  }
  
  tiles4 <- tiles
  tiles4$score <- apply(data.frame(mcols(tiles2)),1,mean)
  rtracklayer::export(tiles4,con=paste(groups[i],".50kb.mean.bedgraph"))
}
