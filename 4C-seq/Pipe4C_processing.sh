# REFERENCE FILE CREATION AND PIPE4C RUN
#
# This script:
#   1) Downloads GRCh38.p12 and prepares an hg38+HPV31 reference FASTA
#   2) Prepares files needed for the BSgenome hg38+HPV31 package
#   3) (BSgenome build and bowtie2 index commands are run in the middle section)
#   4) Runs pipe4C on vp_ecoRI_R1.txt and vp_ecoRI_R2.txt using raw 4C FASTQs
#
# Requirements:
#   - wget
#   - samtools
#   - python
#   - R + pipe4C (Bioconductor)
#   - bowtie2

##############

# REFERENCE FILE CREATION FOR PIPE4C / BSGENOME

cd references
mkdir hg38_hpv31_bsgenome

module load samtools/1.17
module load python

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.p12.genome.fa.gz
gzip -d GRCh38.p12.genome.fa.gz
samtools faidx GRCh38.p12.genome.fa
cut -f1,2 GRCh38.p12.genome.fa.fai > chrom.sizes
python splitRef.py GRCh38.p12.genome.fa
mv GRCh38.p12.genome.fa ..

cp ../HPV31REF_PP706107.fa chrHPV31.fa
# Ensure the FASTA header of chrHPV31.fa is "chrHPV31" (required by hg38_hpv31_seed.txt).
# For example, you can run (uncomment if desired):
# sed -i '1s/.*/>chrHPV31/' chrHPV31.fa

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

# Folder containing the raw 4C FASTQs referenced in vp_ecoRI_R1.txt and vp_ecoRI_R2.txt
DIRFQ=/path/to/4C_fastqs   # EDIT: set to your 4C FASTQ directory

VPFILE=vp_ecoRI_R1.txt
OUTDIR=pipe4c_ecoRI_R1

# NOTE: pipe4C/pipe4C.R is part of the external pipe4C package (Bioconductor).
#       Install pipe4C and update this path if needed.
Rscript pipe4C/pipe4C.R -o $OUTDIR --vpFile=$VPFILE --fqFolder=$DIRFQ

VPFILE=vp_ecoRI_R2.txt
OUTDIR=pipe4c_ecoRI_R2   # EDIT: change if you want a different output directory

Rscript pipe4C/pipe4C.R -o $OUTDIR --vpFile=$VPFILE --fqFolder=$DIRFQ