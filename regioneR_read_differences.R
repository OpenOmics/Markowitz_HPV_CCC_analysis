####################
# Load required packages

if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("GenomicRanges")
}

if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("rtracklayer")
}

if (!requireNamespace("regioneR", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("regioneR")
}


####################
# Functions for prepping files for bootstrapping

expandRegions <- function(grObject, finalSize) {
  # This function will resize all rows in a GenomicRanges object around their
  # midpoints to reach a desired final size.
  library(GenomicRanges)
  midpoint <- floor(width(grObject) / 2)
  start(grObject) <- start(grObject) + midpoint
  end(grObject) <- start(grObject)
  grObject <- shift(grObject, shift=-finalSize/2)
  grObject <- resize(grObject, width=finalSize)
  return(grObject)
}

prepFile <- function(inFile, addFlanks=0) {
  # This function will import a bed or gff file and prepare it for downstream
  # analyses. Steps include importing the data as a GenomicRanges object,
  # sorting it, extending the regions as requested, and collapsing overlapping
  # rows.
  library(GenomicRanges)
  fileData <- rtracklayer::import(inFile)
  fileData <- sort(fileData)
  if (addFlanks != 0) {
    fileData <- regioneR::extendRegions(fileData, extend.start=addFlanks, extend.end=addFlanks)
  }
  return(fileData)
}

getMeansRand <- function(chrsFile, regions, bedgraphFiles, blacklist=NA) {
  chrs <- read.table(chrsFile)
  bedgraphs <- lapply(bedgraphFiles, rtracklayer::import)
  regions$source <- "Original"
  
  if (is.na(blacklist)) {
    rRegions <- randomizeRegions(A=regions, genome=chrs,
                                 allow.overlaps = FALSE)
  } else{
    rRegions <- randomizeRegions(A=regions, genome=chrs, mask=blacklist,
                                 allow.overlaps = FALSE)
  }
  rRegions$source <- "Randomized"
  
  for (i in 1:length(bedgraphFiles)) {
    tmpReal <- sapply(split(regions), meanInRegions, x=bedgraphs[[i]])
    mcols(regions)[i+1] <- tmpReal
    names(mcols(regions))[i+1] <- bedgraphFiles[i]
    tmpRand <- sapply(split(rRegions), meanInRegions, x=bedgraphs[[i]])
    mcols(rRegions)[i+1] <- tmpRand
    names(mcols(rRegions))[i+1] <- bedgraphFiles[i]
  }
  
  outData <- c(regions, rRegions)
  return(as.data.frame(outData))
}

getStats <- function(outData) {
  # wants the direct output from getMeansRand
  IdxOrig <- which(outData$source == "Original")
  IdxRand <- setdiff(1:nrow(outData), IdxOrig)
  for (i in 7:ncol(outData)) {
    Orig <- outData[IdxOrig,i]
    Rand <- outData[IdxRand,i]
    print(names(outData)[i])
    print(wilcox.test(Orig,Rand))
  }
}

###################
## EXAMPLE 1
# This works to get stats where you have the mean "reads" for all the regions of 
# interest and compare it to the mean of random regions of equal size over and over

library(regioneR)

chrs <- read.table("../hg38.fa.sizes")
bedgraph <- rtracklayer::import("9E.50kb.mean.bedgraph")
regions <- prepFile("../Circos/hg38_centromeres.bed")
blacklist <- prepFile("hg38-blacklist.v2.bed")

pt <- permTest(A=regions, ntimes=1000, genome=chrs, mask=blacklist,
               randomize.function=randomizeRegions, allow.overlaps=FALSE,
               evaluate.function=meanInRegions, x=bedgraph)
plot(pt)

# To see actual mean value
pt$meanInRegions$observed

# To get all the permuted mean values
pt$meanInRegions$permuted

###################
## EXAMPLE 2
# This method involves picking a random set of regions once and calculating
# the values of the "reads" for each region separately, producing an output
# table to save for plotting afterward, and a Wilcoxon rank-sum test for stats.
# Note: To make sure you use the same random regions for all sets
#      all bedgraphs must be analyzed at once.

# Regions are the features you want to bootstrap
# Split output file on "source" to for graphics and statistics

library(regioneR)

chrsFile <- "../hg38.fa.sizes"
regions <- prepFile("../Circos/hg38_centromeres.bed")
bedgraphFile1 <- "9E_viewpoint1.50kb.mean.bedgraph"
bedgraphFile2 <- "9E_viewpoint2.50kb.mean.bedgraph"

a <- getMeansRand(chrsFile, regions, 
                  bedgraphFiles=c(bedgraphFile1, bedgraphFile2))

# may give warning messages, if so you can ignore them
# must be done prior to filtering extra columns to work as expected
getStats(a)

# To filter extra columns, optional
a <- a[,c(1:3,6:ncol(a))] 

write.table(x=a,
            file="outfile.txt",
            sep="\t",row.names=F,col.names=T,quote=F)
