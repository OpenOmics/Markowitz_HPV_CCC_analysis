# Markowitz_HPV_CCC_analysis

Citation  
Warburton, Alix et al. “Human papillomavirus genomes associate with active host
chromatin during persistent viral infection.” *PLoS Pathogens* 21(9): e1013454,
2025. doi: 10.1371/journal.ppat.1013454

This repository contains the computational workflows used to analyze HPV31–host
chromatin interactions in CIN612 cells for this paper. It is organized around
two complementary datasets:

- **Hi-C** – processing and analysis of CIN612-9E Hi-C libraries  
- **4C-seq** – processing and analysis of HPV31- and host-centered viewpoints
  in CIN612-6E and CIN612-9E

The code is arranged into:

- `HiC/` – Hi-C processing (reference construction, adapter trimming, HiC-Pro),
  trans-interaction calling (MHiC), and A/B compartment analysis (HiTC).
- `4C-seq/` – 4C-seq reference construction (hg38+HPV31 BSgenome and bowtie2
  index), pipe4C runs, 50 kb binning/normalization, replicate averaging, and
  peak calling.
- Top-level R scripts (e.g. `FigureCreation.R`,
  `regioneR_read_differences.R`) – figure-oriented code and bootstrapping
  analyses that combine processed Hi-C and 4C-seq outputs with external
  epigenomic datasets to generate the panels and enrichment tests shown in the
  manuscript.

The goal of this README is to document **which tools were used**, **how the
files and scripts fit together**, and **how to reproduce the main analyses and
figures at a conceptual level**. It does not attempt to be a full pipeline
wrapper or to reproduce every parameter used on the original cluster.

---

## 1. Data

### Hi-C

- **SM159** – GSM8896999 (CIN612-9E Hi-C)  
- **SM163** – GSM8897000 (CIN612-9E Hi-C)  

In the manuscript SM159 and SM163 are referred to as Replicate 1 and Replicate 2, respectively.
Originally, the samples were named SMC159 and SMC163. This nomencalture may still exist in this repo.

### 4C-seq

- GSM8896992, GSM8896993, GSM8896994 – CIN612-6E 4C-seq (HPV31 viewpoint 3)  
- GSM8896995, GSM8896996, GSM8896997, GSM8896998 – CIN612-9E 4C-seq (HPV31 viewpoint 3)  

These GEO accessions correspond to the 4C-seq datasets used in the manuscript
to map HPV31–host chromatin interactions in CIN612-6E (integrated HPV31) and
CIN612-9E (episomal HPV31) cells.

---

## 2. Software (high-level)

### Core tools

- HiC-Pro (alignment, filtering, contact maps, ICE normalization)  
- pipe4C (4C-seq mapping, fragment counting, and basic QC)  
- bowtie2, samtools, cutadapt, pigz (standard NGS utilities)  
- Python (pandas, matplotlib, seaborn)  
- R + Bioconductor for downstream statistics and figure generation  
  (e.g. MHiC, HiTC, GenomicRanges, BSgenome.Hsapiens.UCSC.hg38, karyoploteR,
  circlize, EnrichedHeatmap, regioneR, rtracklayer, ggplot2)

Exact versions used are noted in comments inside the R and shell scripts.

---

## 3. Key files and what they do

### Hi-C scripts & configs

- `HiC/HiC_processing.sh`  
  - Shows how the combined **hg38 + HPV31** reference was created, adapters
    were trimmed with cutadapt, HiC-Pro was run, and **trans** valid pairs were
    extracted for SM159 and SM163.

- `HiC/TruSeq_and_nextera_adapters.consolidated.fa`  
  - Adapter FASTA used by `HiC_processing.sh` for trimming raw Hi-C reads
    with cutadapt prior to running HiC-Pro.

- `HiC/hicpro_config.txt`  
  - HiC-Pro configuration file (bin sizes, digestion model, ligation site,
    filtering options, output formats). It is used as-is by HiC-Pro and should
    not be edited lightly.

- `HiC/HPV31REF_J04353.fa`  
  - HPV31 FASTA used for the Hi-C portion of the analysis; concatenated to
    GRCh38 to build the combined reference in `HiC_processing.sh`.

- `HiC/HiC_eo.py`  
  - Computes chromosome-association preference (expected/observed
    inter-chromosomal contacts) and produces the heatmaps used for **Figures
    1B–1C** from `SM159.transValidPairs` (and analogously for SM163).
    
- `HiC/HiC_transInteractions_ABcalls.R`  
  - R helpers to run MHiC on HiC-Pro matrices, focus on HPV31–host trans
    interactions, and compute AB (PC1) compartment calls at 500 kb / 1 Mb.

### 4C-seq scripts & configs

- `4C-seq/Pipe4C_processing.sh`  
  - Prepares the combined hg38+HPV31 reference needed for pipe4C (using
    `hg38_hpv31_seed.txt`, `HPV31REF_PP706107.fa`, and `splitRef.py`), builds
    the custom `BSgenome.HsapiensHPV31.custom.hg38` package and bowtie2 index,
    and runs the pipe4C R pipeline on raw 4C-seq FASTQs using
    `vp_ecoRI_R1.txt`, `vp_ecoRI_R2.txt`, and `pipe4c_conf.yml`.

- `4C-seq/hg38_hpv31_seed.txt`  
  - Seed file for `BSgenomeForge` describing the combined hg38+HPV31 genome
    used to build `BSgenome.HsapiensHPV31.custom.hg38`.

- `4C-seq/HPV31REF_PP706107.fa`  
  - HPV31 reference FASTA used in the 4C-seq analyses; merged with hg38 to
    create the combined hg38+HPV31 reference.

- `4C-seq/splitRef.py`  
  - Helper script used by `Pipe4C_processing.sh` to split the hg38 FASTA into
    per-chromosome files and generate `chrom.sizes` for the BSgenome build.

- `4C-seq/pipe4c_conf.yml`  
  - Configuration file for the pipe4C R package, specifying default trimming
    and mapping parameters, restriction enzymes, and the mapping between
    genome identifiers (e.g. `hg38hpv31`) and the corresponding BSgenome /
    bowtie2 entries.

- `4C-seq/vp_ecoRI_R1.txt`, `4C-seq/vp_ecoRI_R2.txt`  
  - Tab-delimited viewpoint definition tables used by pipe4C. Each row
    describes a 4C-seq viewpoint (sample name, restriction enzymes, genome
    identifier, viewpoint chromosome/position, analysis label, and FASTQ file
    name).

- `4C-seq/50kb_bin_normalization_bdg.R`  
  - R script that reads the per-sample pipe4C `.rds` outputs, tiles the genome
    into 50 kb bins, normalizes counts, and writes per-sample 50 kb bedGraph
    files.

- `4C-seq/Average_Viewpoint_50kb.R`  
  - R script that groups 50 kb bedGraphs by condition/viewpoint and writes
    across-replicate mean tracks (e.g. `6E_viewpoint3.50kb.mean.bedgraph`,
    `9E_viewpoint3.50kb.mean.bedgraph`).

- `4C-seq/hist_bedgraphs_50kb_vp3_indSamples_w_6E.R`  
  - R script that performs Poisson-based peak calling on the 50 kb bedGraphs,
    outputting per-sample BED files of significant bins and corresponding
    p-value tracks.

---

### Other

- `FigureCreation.R`  
  - Figure-oriented R script (top-level). Organized into labeled blocks for
    each figure or panel (e.g. Hi-C karyoplots, circos plots, TSS-centered
    profiles, and 4C-seq summary plots). Each block documents the input files
    it expects and which Hi-C / 4C outputs it uses.

---

## 4. Hi-C analysis flow

This is the **high-level order of operations**; exact commands live in the
scripts listed above.

1. **Reference setup**  
   - Download GRCh38.p12 (GENCODE release 28).  
   - Concatenate with `HPV31REF_J04353.fa` to make a combined genome FASTA.  
   - Build chromosome sizes, restriction digest BED, and bowtie2 index, as
     exemplified in `HiC_processing.sh`.  
   - Ensure `hicpro_config.txt` points to the combined reference, sizes and
     digest files.

2. **Pre-processing**  
   - Trim adapters from raw FASTQs for SM159 and SM163 using cutadapt (see the
     example commands and adapter FASTA path in `HiC_processing.sh`).

3. **HiC-Pro**  
   - Run HiC-Pro with `hicpro_config.txt` to generate `.allValidPairs` and ICE-
     normalized contact matrices at the specified resolutions.   

4. **Trans-only valid pairs**  
   - From `SM159.allValidPairs` and `SM163.allValidPairs`, extract inter-
     chromosomal pairs to create `SM159.transValidPairs` and
     `SM163.transValidPairs` (last section of `HiC_processing.sh`). 

5. **Chromosome-association heatmaps (Hi-C)**  
   - Run `HiC_eo.py` to read `SM159.transValidPairs`, compute expected/observed
     inter-chromosomal interaction frequencies, and generate the heatmap PDF
     for Figure 1B (and similarly for SM163 / Figure 1C). 

6. **HPV-host trans calls and compartments**  
   - Use `HiC_transInteractions_ABcalls.R` to:  
     - subset HPV-focused contacts,  
     - run MHiC to call significant HPV-host trans interactions,  
     - compute AB compartment scores from HiC-Pro matrices.

7. **Figures (Hi-C + 4C-seq + epigenomic data)**  
   The script `FigureCreation.R` is organized into labeled blocks 
   by figure/panel. For the Hi-C portion of the paper, the main blocks are:

   - **Figure 1D** – KaryoploteR view of significant HPV31–host trans
     interactions (uses MHiC HPV-significant calls and AB compartment tracks).
   - **Figure 2C** – 1 Mb Pearson correlation heatmap across the genome
     (uses HiC-Pro 1 Mb matrices and derived correlation matrices).
   - **Figure 7D** – Circos plot combining HPV31 4C
     consensus peaks, significant HPV31 trans-interactions from Hi-C, 
     CESC integration hotspots, and W12 super-enhancers to summarize HPV31–host
     connections in CIN612-9E cells.

---

## 5. 4C-seq analysis flow

This is the **high-level order of operations**; exact commands live in the
scripts listed above.

1. **Reference + pipe4C setup**  
   Use `4C-seq/Pipe4C_processing.sh` to prepare the combined hg38+HPV31
   reference required for pipe4C:
   - download GRCh38.p12 and combine it with `HPV31REF_PP706107.fa`,  
   - generate per-chromosome FASTA and `chrom.sizes` using `splitRef.py` and
     `hg38_hpv31_seed.txt`,  
   - build and install `BSgenome.HsapiensHPV31.custom.hg38`,  
   - build a bowtie2 index for the combined genome.

2. **pipe4C: mapping and fragment counting**  
   - Run the pipe4C R pipeline from `Pipe4C_processing.sh` using `vp_ecoRI_R1.txt`,
     `vp_ecoRI_R2.txt`, and `pipe4c_conf.yml` to generate per-sample `.rds`
     files containing fragment-level 4C counts for each viewpoint.

3. **50 kb binning and normalization**  
   - Run `4C-seq/50kb_bin_normalization_bdg.R` on the pipe4C `.rds` outputs to
     tile the genome into 50 kb bins, normalize read counts, and write per-sample
     50 kb bedGraph tracks (`<sample>.50kb.bedGraph`).

4. **Replicate averaging per viewpoint**  
   - Run `4C-seq/Average_Viewpoint_50kb.R` on the 50 kb bedGraphs to group
     samples by condition and viewpoint (e.g. CIN612-6E vs CIN612-9E, HPV31
     viewpoint 3) and write across-replicate mean tracks (e.g.
     `6E_viewpoint3.50kb.mean.bedgraph`, `9E_viewpoint3.50kb.mean.bedgraph`).

5. **Peak calling on 4C profiles**  
   - Run `4C-seq/hist_bedgraphs_50kb_vp3_indSamples_w_6E.R` on the 50 kb
     bedGraphs to perform Poisson-based peak calling and write per-sample peak
     tracks (`<sample>.tem.50kb.bed` and `<sample>.tem.50kb.ppois.bedgraph`).

6. **Figures (Hi-C + 4C-seq + epigenomic data)**  
   The script `FigureCreation.R` is organized into labeled blocks 
   by figure/panel. For the 4C-seq portion of the paper, the main blocks are:

   - **Figure 3** – Circos plot of CIN612-6E HPV31 viewpoint 3
     4C-seq signal across the genome.
   - **Figure 4C** – Circos plot of CIN612-9E
     HPV31 viewpoint 3 4C-seq signal with key genes labeled.
   - **Figure 7D** – Circos plot combining HPV31 4C
     consensus peaks, significant HPV31 trans-interactions from Hi-C, 
     CESC integration hotspots, and W12 super-enhancers to summarize HPV31–host
     connections in CIN612-9E cells.


---
