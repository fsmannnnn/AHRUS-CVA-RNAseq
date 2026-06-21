# AHRUS–CVA bulk RNA-seq & integrative analysis

Code accompanying the manuscript:

> **Cerebrovascular accident as a major AHRUS comorbidity correlates with gut
> dysbiosis-associated CD8+ T cell infiltration.**

This repository contains the **downstream analysis code** for the mouse rectal
bulk RNA-seq (NC / Sham / MCAO, n = 3 per group), the immune-deconvolution and
enrichment analyses, the chemokine visualisation, and the processing of the
public single-cell dataset **GSE264408** used to localise *Ccl5*/*Ccl20* and to
estimate ILC proportions.

---

## 1. Data availability

| Data | Accession / source |
| --- | --- |
| Mouse rectal bulk RNA-seq (NC/Sham/MCAO) |  |
| Mouse fecal 16S rRNA amplicon |  |
| Public single-cell reference (colitis) | GSE264408 |
| Mouse reference genome | Ensembl GRCm38 |

---

## 2. Upstream workflow (performed by the sequencing provider)

### 2.1 Bulk RNA-seq
1. **Adapter & quality trimming — Skewer v0.2.2**
   3′ quality threshold `Q20`, mean-quality filter `Q20`,
   minimum retained length `36bp`.
2. **Quality control — FastQC v0.11.5** (post-trim).
3. **Alignment — STAR v2.5.3** to Ensembl GRCm38 (default parameters unless noted).
4. **Quantification — StringTie v1.3.1**, strand-specific, gene-level counts.

### 2.2 16S rRNA amplicon
1. **Quality filtering — fastp v0.19.6** (sliding window size `4`,
   mean-quality cutoff `Q20`, minimum length `50 bp`).
2. **Read merging — FLASH v1.2.11** (min overlap `10 bp`, max mismatch ratio `0.2`).
3. **Denoising / ASVs — DADA2 within QIIME2 v2020.2.**
4. **Taxonomy — naive Bayes classifier vs. SILVA v138.**
5. Downstream diversity/LEfSe analyses on the Majorbio Cloud Platform.

---

## 3. Downstream workflow (this repository)

All steps are in [`RNAseq_downstream_analysis.R`](RNAseq_downstream_analysis.R),
in the order below.

| # | Step | Key parameters | Figure |
| --- | --- | --- | --- |
| 1 | Read raw counts, collapse duplicate symbols | keep highest mean-expression row | — |
| 2 | Count → **TPM**, low-expression filter | drop genes with **mean TPM ≤ 0.5** | — |
| 3 | Sample selection / relabel | NC1–3, mcao1–3, sham1–3 | — |
| 4 | **PCA** | top-1000 variable genes, log10(TPM+1), scaled | Fig S1A |
| 5 | **DESeq2** DE (MCAO vs Sham; Sham vs NC) | on raw counts; `filterByExpr`; **DEG: p < 0.05 & \|log2FC\| > 1.5** | — |
| 6 | Volcano plot (EnhancedVolcano) | pCutoff 0.05, FCcutoff 1.5 | — |
| 7 | Three-group DEG heatmap | union of DEGs (\|log2FC\|≥1), row Z-score | Fig 4A |
| 8 | Chemokine / inflammation heatmap | curated gene panel, row Z-score | Fig 4J |
| 9 | **CIBERSORT** + **xCell** (IOBR) | CIBERSORT `perm = 500, arrays = FALSE`; xCell `arrays = TRUE` | Fig 4B–C, S1B |
| 10 | **GO / KEGG** ORA (up/down) | clusterProfiler; org.Mm.eg.db / "mmu"; **enrichment input \|log2FC\| > 0.5**; p & q < 0.05 | Fig 4 |
| 11 | **GSEA** | msigdbr (Mus musculus) H / C5-BP / C2-KEGG; rank by log2FC; BH; fgsea; minGSSize 10 (H:20) / maxGSSize 500 (H:1000) | Fig 4D–I |
| 12 | scRNA GSE264408 (Seurat) | LogNormalize; 3000 HVGs; PCs 1:30; UMAP | Fig 4K, S1 I–J |
| 13 | **CIBERSORTx** ILC deconvolution | S-mode, 100 perms (web server) | Fig S1 C–H |

---


## 4. Environment

- **R ≥ 4.4** (Bioconductor 3.19).
- Install from CRAN/Bioconductor; IOBR, COSG and TMEclassifier-style tools are
  installed via `BiocManager`/`devtools` (see script header).
- Random seed fixed (`set.seed(1234)`) for PCA/UMAP reproducibility.

### Quick start
```r
# from the project root
source("RNAseq_downstream_analysis.R")
```

---

## 5. Citing the software

All packages and their developers are listed with versions and DOIs in
[`SOFTWARE_CITATIONS.md`](SOFTWARE_CITATIONS.md).
