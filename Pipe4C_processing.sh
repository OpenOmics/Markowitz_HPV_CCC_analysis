# Prepare files for BSGENOME:
cd references
mkdir hg38_hpv31_bsgenome
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.p12.genome.fa.gz
gzip -d GRCh38.p12.genome.fa.gz
samtools faidx GRCh38.p12.genome.fa
cut -f1,2 GRCh38.p12.genome.fa.fai > chrom.sizes
python splitRef.py GRCh38.p12.genome.fa
mv GRCh38.p12.genome.fa ..

cp ../HPV31REF_PP706107.fa chrHPV31.fa
### manually changed the header to chrHPV31.fa
cd ..

##############

# Build BSGENOME package

module load R/4.4.0
R
library(BSgenome)
BiocManager::install("BSgenomeForge")
forgeBSgenomeDataPkg("hg38_hpv31_seed.txt",replace=T)
q()

R CMD build BSgenome.HsapiensHPV31.custom.hg38
R CMD check BSgenome.HsapiensHPV31.custom.hg38
R CMD INSTALL BSgenome.HsapiensHPV31.custom.hg38

# To install in R as a package:
install.packages("BSgenome.HsapiensHPV31.custom.hg38_0.0.1.tar.gz", repos = NULL, type = "source")

##############

# Build bowtie2 references

mv GRCh38.p12.genome.fa hg38hpv31.fa
cat hg38_hpv31_bsgenome/chrHPV31.fa >> hg38hpv31.fa

module load bowtie/2-2.5.1
bowtie2-build hg38hpv31.fa hg38hpv31
mv hg38hpv31*bt2 bowtie2_ref

###############

# Run pipe4C

module load R/4.4.0
module load bowtie/2-2.5.1
module load samtools/1.17

VPFILE=vp_ecoRI_R1.txt
OUTDIR=pipe4c_ecoRI_R1

Rscript pipe4C/pipe4C.R -o $OUTDIR --vpFile=$VPFILE --fqFolder=$DIRFQ

VPFILE=vp_ecoRI_R2.txt
OUTDIR=/data/NCBR/projects/NCBR-364/pipe4c_ecoRI_R2

Rscript pipe4C/pipe4C.R -o $OUTDIR --vpFile=$VPFILE --fqFolder=$DIRFQ
