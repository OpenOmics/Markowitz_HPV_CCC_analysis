# REFERENCE FILE CREATION FOR HIC-PRO

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.p12.genome.fa.gz
zcat GRCh38.p12.genome.fa.gz > hg38_basic_HPV31.genome.fa
cat HPV31REF_J04353.fa >> hg38_basic_HPV31.genome.fa

module load samtools
samtools faidx hg38_basic_HPV31.genome.fa
cut -f1,2 hg38_basic_HPV31.genome.fa.fai > hg38_basic_HPV31.fa.sizes

module load hicpro/2.11.1
python /hicpro/2.11.1/bin/utils/digest_genome.py -r "mboi" \
          -o hg38_hpv31_mboi.bed hg38_basic_HPV31.genome.fa

module load bowtie/2-2.3.5
bowtie2-build hg38_basic_HPV31.genome.fa hg38_HPV31

#####################

# ADAPTER TRMMING

# NOTE: FASTQ files were too large to run all at once so files were split into 
# 115 chunks for trimming and running through HiC-Pro.

module load cutadapt/1.18

SAMPLE="SM159" # GSM8896999
cutadapt --pair-filter=any --nextseq-trim=2 --trim-n -n 5 -O 5 -q 10,10 \
           -m 35:35 -b file:TruSeq_and_nextera_adapters.consolidated.fa \
           -B file:TruSeq_and_nextera_adapters.consolidated.fa -j 32 \
           -o $SAMPLE.R1.trim.fastq \
           -p $SAMPLE.R2.trim.fastq $SAMPLE.R1.fastq.gz $SAMPLE.R2.fastq.gz 

pigz -p 32 $SAMPLE.R1.trim.fastq
pigz -p 32 $SAMPLE.R2.trim.fastq

SAMPLE="SM163" # GSM8897000
cutadapt --pair-filter=any --nextseq-trim=2 --trim-n -n 5 -O 5 -q 10,10 \
           -m 35:35 -b file:TruSeq_and_nextera_adapters.consolidated.fa \
           -B file:TruSeq_and_nextera_adapters.consolidated.fa -j 32 \
           -o $SAMPLE.R1.trim.fastq \
           -p $SAMPLE.R2.trim.fastq $SAMPLE.R1.fastq.gz $SAMPLE.R2.fastq.gz 

pigz -p 32 $SAMPLE.R1.trim.fastq
pigz -p 32 $SAMPLE.R2.trim.fastq

mv *fastq.gz HicPro

#####################

# EXAMPLE HIC-PRO RUN COMMAND

# NOTE: HiC-Pro was chosen as it was designed to be run step-wise and can easily
# accommodate fastqs that have been split and will need to be consolidated. Code below is
# a basic run example.

HiC-Pro -i HicPro -o HicPro_out -c hicpro_config.txt -s ice_norm

######################

# EXTRACTING TRANS VALID PAIRS

awk '$2!=$5' <SM159.allValidPairs > SM159.transValidPairs
awk '$2!=$5' <SM163.allValidPairs > SM163.transValidPairs
