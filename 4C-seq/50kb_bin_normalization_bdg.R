# Usage:
#   Run from the project root (or adjust paths) after pipe4C has generated
#   per-sample RDS files in:
#     pipe4c_ecoRI_R1/RDS/*.rds
#     pipe4c_ecoRI_R2/RDS/*.rds
#
#   Output:
#     One 50 kb bedGraph per unique 4C sample:
#       <sample>.50kb.bedGraph
#     in the current working directory.
#
# Requirements:
#   - R + packages: GenomicRanges, BSgenome.Hsapiens.UCSC.hg38, rtracklayer
#   - The combined hg38+HPV31 BSgenome is NOT required here; only hg38.

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

norm4C <- function( readsGR, nReads=1e6, nTop=2, wSize=21 ) {
  readsGR$normReads <- 0
  sumTop <- sum( -sort( -readsGR$reads )[ 1:nTop ] )
  wNorm <- nReads/( sum( readsGR$reads )-sumTop )
  readsGR$normReads <- wNorm*readsGR$reads
  #readsGR$norm4C <- caTools::runmean( x=readsGR$normReads, k=wSize, endrule="mean" )
  return( readsGR )
}

binReads <- function(readsGR, res) {
  tmp <- seqinfo(Hsapiens)
  seqlevels(tmp) <- seqlevels(readsGR)[1:24]
  tiles <- tileGenome(tmp, tilewidth=res,
                    cut.last.tile.in.chrom=TRUE)
  ov1 <- countOverlaps(tiles, readsGR)
  tiles$nRE <- ov1
  readsGR2 <- readsGR[which(readsGR$reads > 0)]
  ov2 <- findOverlaps(tiles, readsGR2)
  sums <- vector(length=length(tiles))
  for (i in 1:length(tiles)) {
       tmp <- which(queryHits(ov2) == i)
       sums[i] <- sum(readsGR2$normReads[subjectHits(ov2)[tmp]])
  }
  tiles$sumNormReads <- sums
  tiles <- tiles[which(tiles$sumNormReads > 0)]
}


R1files <- list.files(path="pipe4c_ecoRI_R1/RDS",full.names=T)
R2files <- list.files(path="pipe4c_ecoRI_R2/RDS",full.names=T)
R1samples <- gsub("pipe4c_ecoRI_R1/RDS/","",gsub(".rds","",R1files))
R2samples <- gsub("pipe4c_ecoRI_R2/RDS/","",gsub(".rds","",R2files))

uniqueSamples <- sort(unique(c(R1samples,R2samples)))
allFiles <- sort(c(R1files, R2files))

outData2 <- vector("list",length=length(uniqueSamples))                

for (i in 1:length(uniqueSamples)) {
files <- grep(uniqueSamples[i], allFiles, value=T)

if (length(files) == 1) {
  inData <- readRDS(files)
  outData <- inData$reads
  mcols(outData) <- NULL
  mcols(outData)$reads <- inData$reads$reads
} else {
  inData <- lapply(files, readRDS)
  outData <- inData[[1]]$reads
  mcols(outData) <- NULL
  mcols(outData)$reads <- inData[[1]]$reads$reads + inData[[2]]$reads$reads
}
 outData2[[i]] <- norm4C(outData)
}

for (i in 1:length(uniqueSamples)) {
  print(i)
  a <- binReads(outData2[[i]],5e4)
  a$score <- a$sumNormReads
  
  outFile <- paste0(uniqueSamples[i], ".50kb.bedGraph")
  rtracklayer::export(a, outFile)
}