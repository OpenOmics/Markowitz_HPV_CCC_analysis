#################
# Figure 1D: KaryoploteR of HiC trans interaction calls

library(GenomicRanges)
library(karyoploteR)

# Bed file used here was all peaks from SM159_1000000_noYM_mhic_output_table_HPVsig.txt
# and SM163_1000000_noYM_mhic_output_table_HPVsig.txt
# Sources:
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8896999
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8897000

a <- read.table("sigHPV_noYM_1Mb.bed")
names(a) <- c("chr","start","end","name","score")
Mb <- makeGRangesFromDataFrame(a, keep.extra.columns=T)
Mb163 <- Mb[grep("163",Mb$name)]
Mb159 <- Mb[grep("159",Mb$name)]

kp <- plotKaryotype(genome="hg38", plot.type=1)
kpPlotRegions(kp,reduce(Mb163),r0=0,r1=0.2,col="royalblue4")
kpPlotRegions(kp,reduce(Mb159),r0=0.3,r1=0.5,col="darkorange")

#################
# Figure 2C: 1Mb pearson correlation heatmap

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)

allFiles <- list.files(pattern="bedgraph")
allData <- lapply(allFiles,rtracklayer::import)

tmp <- seqinfo(Hsapiens)
seqlevels(tmp) <- seqlevels(allData[[1]])[1:24]
tiles <- tileGenome(tmp, tilewidth=50000,
                     cut.last.tile.in.chrom=TRUE)

ovs <- lapply(allData,function(x) {findOverlaps(x,tiles)})


M1 <- data.frame(data.frame(tiles)[,1:3], matrix(ncol=length(allFiles),nrow=length(tiles),data=0))

for (i in 1:length(allFiles)) {
   M1[subjectHits(ovs[[i]]), (i+3)] <- allData[[i]]$score
}

names(M1)[4:ncol(M1)] <- gsub(".mean.bedgraph","",allFiles)
names(M1)[4:ncol(M1)] <- c("VP1","VP2","VP3")
M1corr <- cor(M1[1:(nrow(M1)-1),4:ncol(M1)],method="pearson")
M1corrM <- reshape2::melt(M1corr)

ggplot(data = M1corrM, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  geom_text(aes(label = round(value, 2))) +
  scale_y_discrete(limits=rev) +
  scale_fill_distiller(palette = "YlOrRd",limits=c(0,1),direction = 1) +
  theme_bw() 


#################
# Figure 3: 6E circos

library(circlize)
library(GenomicRanges)

# Sources
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8896992
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8896993
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8896994

bdg6E <- "6E_viewpoint3.50kb.mean.bedgraph"
bdg <- rtracklayer::import(bdg6E)
seqlevels(bdg) <- gsub("chr","",seqlevels(bdg))
bdg2 <- data.frame(bdg)[,c(1,2,3,6)]
  

  circos.genomicInitialize(chrLen2, plotType = NULL)
  
  circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.05,
                         panel.fun = function(x, y) {
                           xlim = get.cell.meta.data("xlim")
                           chr = get.cell.meta.data("sector.index")
                           circos.text(mean(xlim), 0.5, labels = chr, facing = "clockwise", niceFacing = TRUE)
                         })
  
  circos.trackPlotRegion(ylim = c(0, 1), bg.border = col, bg.col = col, panel.fun = function(x, y) {
    chr = get.cell.meta.data("sector.index")
    circos.axis(h = "bottom", labels = NULL, sector.index = chr, direction = "inside")
  }, track.height = 0.05)
  
  circos.genomicTrack(bdg2,ylim=c(0,5000),
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region, value, area = TRUE, col = NA, border = "black", ...)
                      }, bg.border = NA, track.height = 0.3)

#################
# Figure 4C: 9E circos with gene labels

library(circlize)
library(GenomicRanges)
library(regioneR) # to convert data.frame to gRanges

# Sources
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8896995 
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8896996
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8896997
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8896998

chrLen <- read.table("hg38.HPV.sizes",sep="\t")
chrLen2 <- data.frame(chr=chrLen$V1, start=1, end=chrLen$V2)
chrLen2$chr <- gsub("chr","",chrLen2$chr)

col <- paste0("grey",seq(5,by=4,length.out=24))

# coordinates grabbed from Ensembl GRCh38.p14
geneList <- data.frame(chr=c(8,1,2),
                       start=c(127735434, 23506438, 201233443),
                       end=c(127742951, 23531233, 201361836),
                       name=c("MYC","E2F2","CASP8"))
geneList <- toGRanges(geneList)

# For 9E vp3
bdg9E <- "9E_viewpoint3.50kb.mean.bedgraph"
bdg <- rtracklayer::import(bdg9E)
seqlevels(bdg) <- gsub("chr","",seqlevels(bdg))
bdg2 <- data.frame(bdg)[,c(1,2,3,6)]

circos.genomicInitialize(chrLen2, plotType = NULL)

circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.05,
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         chr = get.cell.meta.data("sector.index")
                         circos.text(mean(xlim), 0.5, labels = chr, facing = "clockwise", niceFacing = TRUE)
                       })

circos.trackPlotRegion(ylim = c(0, 1), bg.border = col, bg.col = col, panel.fun = function(x, y) {
  chr = get.cell.meta.data("sector.index")
  circos.axis(h = "bottom", labels = NULL, sector.index = chr, direction = "inside")
}, track.height = 0.05)

circos.genomicTrack(bdg2,ylim=c(0,500),
  panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, area = TRUE, col = NA, border = "black", ...)
}, bg.border = NA, track.height = 0.3)

circos.genomicLabels(geneList, labels.column = 4, side = "inside")

#################
# Figure 7D: 9E circos with linkages

library(circlize)
library(GenomicRanges)

# Sources:
#    S1 Table. Significant HPV31 trans-interactions mapped in CIN612-9E cells by HiC (mhic_sigHPVtrans_hg38_FINAL.bed)
#    S6 Table. HPV31 4C consensus peaks defined in CIN612-9E cells (9Evp3_tem_50kb_inMost.FIXED.bed)
#    S8 Table. Super-enhancers profiled in W12 HPV16-positive cervical keratinocytes [Warburton, 2021]
#    S10 Table. Integration hotspots defined in CESC [Warburton, 2021]

CESC <- rtracklayer::import("Final CESC hotspots_hg38_LiftOver.bed")
seqlevels(CESC) <- gsub("chr","",seqlevels(CESC))
CESC <- data.frame(CESC)

Enhancer <- rtracklayer::import("W12 SE_hg38_LiftOver 1.bed")
seqlevels(Enhancer) <- gsub("chr","",seqlevels(Enhancer))
Enhancer <- data.frame(Enhancer)

chrLen <- read.table("hg38.HPV.sizes",sep="\t")
chrLen2 <- data.frame(chr=chrLen$V1, start=1, end=chrLen$V2)
chrLen2$chr <- gsub("chr","",chrLen2$chr)

trans <- rtracklayer::import("9Evp3_tem_50kb_inMost.FIXED.bed")
seqlevels(trans) <- gsub("chr","",seqlevels(trans))

transA <- trans
end(transA) <- end(transA) + 1e6

chrLen3 <- data.frame(chr=rep(chrLen2$chr[24],length(trans)),
                      start=rep(chrLen2$start[24],length(trans)),
                      end=rep(chrLen2$end[24],length(trans)))

HiC <- rtracklayer::import("mhic_sigHPVtrans_hg38_FINAL.bed")
seqlevels(HiC) <- gsub("chr","",seqlevels(HiC))

trans2 <- subsetByOverlaps(trans,HiC)
end(trans2) <- end(trans2) + 1e6

col <- paste0("grey",seq(5,by=4,length.out=24))

circos.genomicInitialize(chrLen2, plotType = NULL)

circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.05,
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         chr = get.cell.meta.data("sector.index")
                         circos.text(mean(xlim), 0.5, labels = chr, facing = "clockwise", niceFacing = TRUE)
                       })

circos.trackPlotRegion(ylim = c(0, 1), bg.border = col, bg.col = col, panel.fun = function(x, y) {
  chr = get.cell.meta.data("sector.index")
  circos.axis(h = "bottom", labels = NULL, sector.index = chr, direction = "inside")
}, track.height = 0.05)

circos.genomicDensity(CESC, col = c("brown"), track.height = 0.1)
circos.genomicDensity(Enhancer, col = c("darkgreen"), track.height = 0.1)

circos.genomicLink(region1=data.frame(transA)[,1:3], region2=chrLen3, 
                   col = "deepskyblue", border=NA)

circos.genomicLink(data.frame(trans2)[,1:3], chrLen3[1:length(trans2),], 
                   col = "red")

#################
# Figure 8A: Profile plot around TSS

# Before running this script, individual bigwig files were averaged using the 
# deeptools function bigwigAverage

library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(EnrichedHeatmap)

# Sources:
#   NHEK DSBCapture BREAK-seq (GSM2068755 and GSM2068756) 
#   20863 H3K27ac ChIP-seq (GSM5550313 and GSM5550314)
#   9E ATAC-seq (GSM8898200 and GSM8898201)

bws <- c("BREAK_n12_avg.bw", "G1.avg.genrich.RPM.bw", "H3K27ac_20863.avg.Q5DD.RPGC.inputnorm.bw")
bwNames <- c("BREAK", "ATAC", "H3K27ac")

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gtf <- import("gencode.v19.annotation.gtf.gz")
gtf <- gtf[which((gtf$type == "gene") & (gtf$gene_type == "protein_coding"))]

# To orient all start codons to the right of the TSS
gtf_starts <- gtf
start(gtf_starts) <- ifelse( strand(gtf) == '+', start(gtf), end(gtf) )
end(gtf_starts) <- start(gtf_starts)

for (i in 1:length(bws)) {
  print(i)
  signal <- import(bws[i])
  tmp <- normalizeToMatrix(signal, gtf_starts,
                               extend=c(5000, 5000), w=25,
                               mean_mode="weighted", value_column="score")
  tmpMean <- colMeans(tmp, na.rm = T)
  tmpMean2 <- data.frame(SampleID=bwNames[i], Position=seq(-199, 200), Mean=tmpMean)
  tmpMean3 <- tmpMean2[,2:3]
  names(tmpMean3)[2] <- bwNames[i]
  if (i == 1) {
     sepMeans2 <- tmpMean3
  } else {
     sepMeans2 <- cbind(sepMeans2, tmpMean3) 
  }
}

# getting rid of the extra "Position" columns
sepMeans2 <- sepMeans2[,c(1,2,4,6)]

colrs <- c("black", "#005AB5", "#E66100")
coeff=2
shiftVal=1.2

ggplot(sepMeans2, aes(x=Position)) +
  geom_line(aes(y=H3K27ac), color=colrs[2],size=1.25) +
  theme_classic() +
   labs(x="") +
   scale_x_continuous(limits=c(-199,200),
                     breaks = c(-199, 0, 200),
                     labels = c("-5kb", "TSS","5kb")) +
   geom_line(aes(y=BREAK/coeff-shiftVal), color=colrs[3], size=1.25) +
   scale_y_continuous(name="ATAC/H3K27ac",
         sec.axis = sec_axis(~.*coeff+shiftVal,name="BREAK")) +
   theme(axis.title.y.right=element_text(color= colrs[3])) +
  geom_line(aes(y=ATAC*3), color=colrs[1],size=1.25)
