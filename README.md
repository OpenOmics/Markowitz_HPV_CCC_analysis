# Markowitz_HPV_CCC_analysis

Citation
--------
Warburton, Alix et al. “Human papillomavirus genomes associate with active host chromatin during persistent viral infection.” PLoS pathogens vol. 21,9 e1013454. 2 Sep. 2025, doi:10.1371/journal.ppat.1013454


## HiC

Minimal run order 
-----------------
Follow these steps in sequence to reproduce Hi‑C trans‑interaction results used in this project:

1. Prepare combined reference
   - Use HiC/HPV31REF_J04353.fa together with a GRCh38 FASTA to create a merged reference for alignment.

2. Create restriction digest and aligner indexes
   - Consult HiC/HiC_processing.sh for the digest and index steps; these produce files used by HiC‑Pro.

3. Preprocess raw reads
   - Use the trimming/preparation steps documented in HiC/HiC_processing.sh (TruSeq_and_nextera_adapters.consolidated.fa referenced for this step).

4. Configure HiC‑Pro
   - See HiC/hicpro_config.txt for chosen HiC-Pro run parameters.

5. Run HiC‑Pro
   - Execute HiC‑Pro (stepwise or single run) to produce per‑sample outputs including .allValidPairs and normalized matrices.

6. Extract inter‑chromosomal (trans) valid pairs
   - Produce per‑sample trans validPairs files from HiC‑Pro outputs (See HiC/HiC_processing.sh for details).

7. Call and summarize trans interactions
   - Run HiC/HiC_transInteractions_ABcalls.R to generate MHic outputs and filtered files used downstream (e.g., *_mhic_output_table_HPVsig.txt).

8. Optional helpers / formatting
   - Use HiC/HiC_eo.py for small text/format conversions if needed.

Software & package versions
---------------------------
System / bioinformatics tools
- HiC‑Pro 2.11.1
- Bowtie2 2.3.5
- cutadapt 1.18
- samtools

R and Bioconductor packages
- MHiC v0.5.0 (R 3.6.1)
- HiTC v1.48.0 (R 4.4.1)
- GenomicFeatures
- AnnotationDbi
- data.table
- GenomicRanges
- GenomeInfoDb
