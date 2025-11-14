# MHiC version 0.5.0 (R/3.6.1)

mhic_run <- function(root, res) {
        library(MHiC)
        
        readsFile = paste0(root,".HPV.matrix")
        digestFile = paste0(root,"_abs.bed")
        sampleName = paste0(root,"_HPV_mhic")
        
        Output <- MHiC(readsFile, digestFile, outdir=sampleName, 
                sample_name=sampleName, tools_name="HiC_PRO", 
                res = as.integer(res), cistrans = "trans", parallel=TRUE, cores=2)
        
        a <- Output[which(Output$qvalue < 0.1),]
        write.table(a, paste0(root,"_mhic_output_table_sig.txt"),
                quote=F, sep="\t", row.names=F)
        
        a <- Output[which((Output$chr1 == "chrHPV31REF") | 
                                        (Output$chr2 == "chrHPV31REF") ),]
        write.table(a, paste0(root,"_mhic_output_table_HPV.txt"),
                quote=F, sep="\t", row.names=F)
        
        b <- a[which(a$qvalue < 0.1),]
        write.table(b, paste0(root,"_mhic_output_table_HPVsig.txt"),
                quote=F, sep="\t", row.names=F)
}

rmYM <- function(root) {
        bedFile <- paste0(root,"_abs.bed")
        outBedFile <- paste0(root,"_noYM_abs.bed")
        matrixFile <- paste0(root,".matrix")
        outmatrixFile <- paste0(root,"_noYM.matrix")
        
        bedData <- read.table(bedFile)
        
        IDs2remove <- which( (bedData[,1] == "chrY") | (bedData[,1] == "chrM") )
        
        newBed <- bedData[which(!(bedData[,4] %in% IDs2remove)),]
        write.table(newBed, outBedFile, sep="\t", quote=F, col.names=F, row.names=F)
        
        matrixData <- read.table(matrixFile)
        newMatrix <- matrixData[which( !( (matrixData[,1] %in% IDs2remove) |  
                                          (matrixData[,2] %in% IDs2remove) ) ), ]
        write.table(newMatrix, outmatrixFile, sep="\t", quote=F, col.names=F, row.names=F)
}

HPVselect <- function(root) {
        bedFile <- paste0(root,"_noYM_abs.bed")
        matrixFile <- paste0(root,"_noYM.matrix")
        outmatrixFile <- paste0(root,"_noYM.HPV.matrix")

        bedData <- read.table(bedFile)

        IDs2keep <- bedData[which(bedData[,1] == "HPV31REF"),4]

        matrixData <- read.table(matrixFile)
        newMatrix <- matrixData[which( (matrixData[,1] %in% IDs2keep) |
                                                                        (matrixData[,2] %in% IDs2keep) ), ]
        write.table(newMatrix, outmatrixFile, sep="\t", quote=F, col.names=F, row.names=F)
}


# To create:
GSM8896999_SM159_1000000_noYM_mhic_output_table_sig.txt.txt.gz	

rmYM("SM163_1000000")
mhic_run("SM163_1000000_noYM",1000000)
rmYM("SM159_1000000")
mhic_run("SM159_1000000_noYM",1000000)
HPVselect("SM159_1000000")
mhic_run("SM159_1000000_noYM",1000000)


###############################
# TXDB creation
# Done by: Brittany Dulek

library(GenomicFeatures)
library(AnnotationDbi)
library(data.table)

chrom_file <- data.table(read.table('hg38.fa.sizes', col.names = c("chromosome", "size"))) # chromosome sizes file

seq.info <- Seqinfo(seqnames = chrom_file$chromosome, 
                    seqlengths = chrom_file$size, 
                    isCircular = c(as.vector(rep(FALSE, 24)), TRUE), # 24 chromosomes not circular.  MT is circular.
                    genome = "GRCh38") # Genome

txdb <- makeTxDbFromGFF("genes_protein_coding.gtf", # name of gtf, only protein coding were used for ATAC-seq
                        format = "gtf", 
                        dataSource = "ATAC-seq Pipeline GTF Gencode 28", # Notate gtf source
                        organism = "Homo sapiens", # organism
                        chrominfo = seq.info)

saveDb(txdb, "ATACseq_pipeline_protein-coding_hg38.txdb") # name of txdb to be created

###############################
# HiTC version 1.48.0 (R/4.4.1)

library(HiTC)
library(GenomicFeatures)

SMC159 <- importC('SM159_500000_iced.matrix', 'SM159_500000_abs.bed')

SMC163 <- importC('SM163_500000_iced.matrix', 'SM163_500000_abs.bed')

txdb <- loadDb("ATACseq_pipeline_protein-coding_hg38.txdb")
geneannot <- genes(txdb)
geneannot <- sort(geneannot)

ToSelf <- paste0(seqlevels(geneannot),seqlevels(geneannot))[1:23]
Idx <- which(names(SMC159) %in% ToSelf)

pc <- vector(mode="list", length=length(Idx))
for (i in 1:length(Idx)) {
  print(names(SMC159)[Idx[i]])
  pc[[i]] <- pca.hic(SMC159[[Idx[i]]], normPerExpected=TRUE, method="loess", npc=1, gene.gr=geneannot)
}

pcB <- lapply(pc, function(x) {x$PC1})
pcC <- do.call("c",pcB)

write.table(data.frame(pcC),"SMC159_ABcalls_per500kb.txt",quote=F,sep="\t")

###

ToSelf <- paste0(seqlevels(geneannot),seqlevels(geneannot))[1:23]
Idx <- which(names(SMC163) %in% ToSelf)

pc <- vector(mode="list", length=length(Idx))
for (i in 1:length(Idx)) {
  print(names(SMC163)[Idx[i]])
  pc[[i]] <- pca.hic(SMC163[[Idx[i]]], normPerExpected=TRUE, method="loess", npc=1, gene.gr=geneannot)
}

pcB <- lapply(pc, function(x) {x$PC1})
pcC <- do.call("c",pcB)

write.table(data.frame(pcC),"SMC163_ABcalls_per500kb.txt",quote=F,sep="\t")




