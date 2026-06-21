# =============================================================================
# Downstream bulk RNA-seq & integrative analysis
# Project : Cerebrovascular accident as a major AHRUS comorbidity correlates
#           with gut dysbiosis-associated CD8+ T cell infiltration
# Tissue  : Mouse rectum, three groups (NC / Sham / MCAO), n = 3 per group
# Scope   : Count->TPM, DESeq2 DE, immune deconvolution (CIBERSORT/xCell),
#           GO/KEGG/GSEA, chemokine heatmap, and public scRNA (GSE264408)
#           processing for CCL5/CCL20 source + CIBERSORTx ILC deconvolution.
#
# NOTE ON UPSTREAM: adapter/quality trimming (Skewer), QC (FastQC), alignment
#   (STAR) and quantification (StringTie) were performed by the sequencing
#   provider. See README.md for the upstream workflow and parameters.
#
# Reproducibility: run sessionInfo() at the end; package versions are pinned
#   in README.md. Replace the placeholder paths under "PATHS" with your own.
# =============================================================================


# ----------------------------- 0. Packages ----------------------------------
suppressPackageStartupMessages({
  # data handling
  library(tidyverse)      # dplyr / tidyr / ggplot2 / stringr / readr
  library(data.table)     # fread / fwrite
  library(readxl); library(openxlsx)
  # differential expression
  library(DESeq2)
  library(edgeR)          # filterByExpr
  library(EnhancedVolcano)
  # dimensionality reduction / visualization
  library(factoextra)
  library(genefilter)
  library(pheatmap)
  library(ComplexHeatmap); library(circlize); library(RColorBrewer)
  # immune deconvolution
  library(IOBR)           # deconvo_tme(): CIBERSORT, xCell
  # enrichment
  library(clusterProfiler); library(org.Mm.eg.db); library(msigdbr)
  library(enrichplot); library(GseaVis)
  # single cell
  library(Seurat)
})

set.seed(1234)

# ----------------------------- PATHS ----------------------------------------
# Expected layout (relative to project root):
#   data/        raw count matrix, gene length table, group file
#   results/     DEG tables, deconvolution outputs
#   results/figures/  publication figures
proj_root <- "."                      # <- set to your project root
dir_data  <- file.path(proj_root, "data")
dir_res   <- file.path(proj_root, "results")
dir_fig   <- file.path(proj_root, "results", "figures")
dir.create(dir_res, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

# group colours used throughout
grp_cols <- c(NC = "#6699CC", sham = "#99CC99", mcao = "#FF9999")


# ======================= 1. Read raw count matrix ===========================
# data.txt: first columns are annotation incl. a 'Symbol' column, then samples.
mrna_counts <- read.delim(file.path(dir_data, "data.txt"),
                          header = TRUE, sep = "\t", quote = "",
                          stringsAsFactors = FALSE, fileEncoding = "UTF-8")
mrna_counts <- as.data.frame(mrna_counts)

# collapse duplicated gene symbols: keep the row with the highest mean expression
idx          <- order(rowMeans(mrna_counts[, -c(1, 2, 3)]), decreasing = TRUE)
expr_ordered <- mrna_counts[idx, ]
mrna_counts  <- expr_ordered[!duplicated(expr_ordered$Symbol), ]
mrna_counts$Symbol <- toupper(mrna_counts$Symbol)
rownames(mrna_counts) <- mrna_counts$Symbol

# sample group table: columns 'sample' and 'group'
mrna_sample <- read.table(file.path(dir_data, "group.txt"),
                          header = TRUE, check.names = FALSE)
rownames(mrna_sample) <- mrna_sample$sample


# ================ 2. Count -> TPM and low-expression filtering ==============
# gene length from Mus_musculus.GRCm38.102.gtf (mm10)
gene_length <- read.csv(file.path(dir_data, "mouse_gene_length.csv"),
                        header = TRUE, sep = ",", stringsAsFactors = FALSE)
gene_length$external_gene_name <- toupper(gene_length$external_gene_name)
colnames(gene_length) <- c("Geneid","Symbol","Chr","Start","End","Strand","Length")
gene_length <- gene_length[, c("Symbol","Length")]

gene_count <- inner_join(gene_length, mrna_counts, by = "Symbol")

# TPM conversion
tpm_value <- gene_count
for (i in 3:ncol(gene_count)) {
  tpm_value[, i] <- round(
    (gene_count[, i] * 1e3 * 1e6) /
      (gene_count[, "Length"] * sum(gene_count[, i] * 1e3 / gene_count[, "Length"])),
    3)
}
tpm_result <- tpm_value[, -2]                    # drop Length column

# remove genes not expressed in any sample, then collapse duplicates again
zero_pct   <- rowMeans(tpm_result[, -1] == 0)
tpm_result <- tpm_result[zero_pct < 1, ]
idx        <- order(rowMeans(tpm_result[, -1]), decreasing = TRUE)
tpm_result <- tpm_result[idx, ]
tpm_result <- tpm_result[!duplicated(tpm_result$Symbol), ]
rownames(tpm_result) <- tpm_result$Symbol
tpm_result <- tpm_result[, -1]

# low-expression filter: keep genes with mean TPM > 0.5
tpm_result <- tpm_result[rowMeans(tpm_result) > 0.5, ]
fwrite(data.frame(gene_name = rownames(tpm_result), tpm_result),
       file.path(dir_data, "data_tpm_filtered.txt"))
cat("Filtered TPM matrix:", nrow(tpm_result), "genes x", ncol(tpm_result), "samples\n")


# =================== 3. Sample selection & grouping =========================
# keep one representative set of NC / MCAO / Sham (n = 3 each)
tpm_filt <- tpm_result %>%
  dplyr::select(starts_with("NC"), mcao2, mcao3, mcao5, starts_with("sham"))

mrna_sample_filt <- mrna_sample[colnames(tpm_filt), ]

new_names <- c("NC1","NC2","NC3","mcao1","mcao2","mcao3","sham1","sham2","sham3")
colnames(tpm_filt)            <- new_names
rownames(mrna_sample_filt)    <- new_names
mrna_sample_filt$sample       <- new_names
stopifnot(all(colnames(tpm_filt) == rownames(mrna_sample_filt)))


# ============================== 4. PCA (Fig S1A) ============================
data_fixed <- as.matrix(tpm_filt); data_fixed[data_fixed < 0] <- 0
rv     <- genefilter::rowVars(data_fixed)
select <- order(rv, decreasing = TRUE)[seq_len(1000)]      # top-1000 variable genes
pca_in <- cbind(t(log10(data_fixed[select, ] + 1)), mrna_sample_filt)
expr_pca <- prcomp(pca_in[, 1:1000], scale. = TRUE, center = TRUE)

PCA <- fviz_pca_ind(expr_pca, label = "none", geom.ind = "point",
                    habillage = mrna_sample_filt$group,
                    addEllipses = TRUE, ellipse.type = "norm",
                    ellipse.level = 0.88, mean.point = FALSE,
                    palette = unname(grp_cols)) +
  theme_bw(base_size = 14) + ggtitle("PCA") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_fixed() + theme(aspect.ratio = 1, legend.title = element_blank())
ggsave(file.path(dir_fig, "PCA_Plot.pdf"), PCA, width = 4, height = 4)
ggsave(file.path(dir_fig, "PCA_Plot.png"), PCA, width = 4, height = 4, dpi = 300)


# =============== 5. DESeq2 differential expression ==========================
# DE is run on RAW COUNTS (not TPM). Build a count matrix matching tpm_filt.
counts_filt <- mrna_counts[rownames(tpm_filt), colnames(mrna_counts) %in% new_names |
                             colnames(mrna_counts) %in% mrna_sample_filt$sample]
# If your raw count columns use the original sample names, re-map them here so
# that columns are exactly: NC1..NC3, mcao1..mcao3, sham1..sham3.

run_deseq2 <- function(counts, meta, contrast) {
  meta$group <- factor(meta$group)
  dds  <- DESeqDataSetFromMatrix(round(as.matrix(counts)), colData = meta, design = ~ group)
  keep <- filterByExpr(DGEList(counts = as.matrix(counts)),
                       model.matrix(~ 0 + meta$group))      # edgeR default expr filter
  dds  <- dds[keep, ]
  dds  <- DESeq(dds, quiet = TRUE)
  res  <- results(dds, contrast = c("group", contrast[1], contrast[2]))
  res  <- as.data.frame(res[order(res$padj), ])
  res  <- na.omit(res); res$gene <- rownames(res)
  res
}

meta_all <- mrna_sample_filt
DEG_mcao_sham <- run_deseq2(counts_filt, meta_all, c("mcao","sham"))
DEG_sham_nc   <- run_deseq2(counts_filt, meta_all, c("sham","NC"))
write.csv(DEG_mcao_sham, file.path(dir_res, "DESeq2_DEG_mcao_sham.csv"), row.names = FALSE)
write.csv(DEG_sham_nc,   file.path(dir_res, "DESeq2_DEG_sham_nc.csv"),   row.names = FALSE)

# DEG definition reported in the manuscript: p < 0.05 & |log2FC| > 1.5
deg_call <- function(df, lfc = 1.5)
  df %>% filter(!is.na(padj), pvalue < 0.05, abs(log2FoldChange) > lfc)


# ---------------------------- 6. Volcano plot -------------------------------
EnhancedVolcano(DEG_mcao_sham,
                lab = DEG_mcao_sham$gene, x = "log2FoldChange", y = "padj",
                title = "MCAO vs. Sham", subtitle = NULL,
                col = c("grey30","#458c45","#3c61a3","#d4372a"),
                pCutoff = 0.05, FCcutoff = 1.5,
                pointSize = 1.0, labSize = 2.5,
                gridlines.major = FALSE, gridlines.minor = FALSE)


# ============= 7. Three-group DEG heatmap (Fig 4A) ==========================
deg_genes <- union(
  (DEG_mcao_sham %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1))$gene,
  (DEG_sham_nc   %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1))$gene)
deg_genes <- intersect(rownames(tpm_filt), deg_genes)

grp_order <- mrna_sample_filt %>% arrange(factor(group, levels = c("NC","sham","mcao")))
mat <- log2(as.matrix(tpm_filt[deg_genes, rownames(grp_order)]) + 1)

anno_col <- data.frame(Group = grp_order$group, row.names = rownames(grp_order))
pheatmap(mat, scale = "row",
         cluster_rows = TRUE, cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = TRUE,
         annotation_col = anno_col,
         annotation_colors = list(Group = grp_cols),
         color = colorRampPalette(c("navy","white","firebrick3"))(100),
         border_color = NA,
         gaps_col = c(sum(grp_order$group == "NC"),
                      sum(grp_order$group %in% c("NC","sham"))),
         filename = file.path(dir_fig, "DEGs_heatmap_three_groups.pdf"),
         width = 8, height = 10)


# ====== 8. Chemokine / inflammation gene heatmap (Fig 4J) ===================
# gene panel read from an external curated list (one column 'Symbol')
infl <- read_excel(file.path(dir_data, "heatmap_genes_filtered.xlsx"))
infl_mat <- tpm_filt[intersect(toupper(infl$Symbol), rownames(tpm_filt)),
                     rownames(grp_order), drop = FALSE]
infl_scaled <- t(scale(t(as.matrix(infl_mat))))

ha <- HeatmapAnnotation(group = grp_order$group,
                        col = list(group = grp_cols),
                        simple_anno_size = unit(0.5, "cm"))
Heatmap(infl_scaled, name = "Z-score",
        col = colorRampPalette(c("#313695","#FFFFFF","#A50026"))(100),
        top_annotation = ha,
        cluster_rows = TRUE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))


# ====== 9. Immune deconvolution: CIBERSORT + xCell (Fig 4B-C, S1B) ==========
decon_in <- as.matrix(tpm_filt)        # genes x samples, gene symbols as rownames

im_cibersort <- deconvo_tme(eset = decon_in, method = "cibersort",
                            arrays = FALSE, perm = 500)
im_xcell     <- deconvo_tme(eset = decon_in, method = "xcell",
                            arrays = TRUE)
write.xlsx(im_cibersort, file.path(dir_res, "im_cibersort_NC_sham_mcao.xlsx"))
write.xlsx(im_xcell,     file.path(dir_res, "im_xcell_NC_sham_mcao.xlsx"))

# CIBERSORT stacked bar (Fig 4B)
res_cibersort <- cell_bar_plot(input    = im_cibersort,
                               features = colnames(im_cibersort)[2:23],
                               title    = "CIBERSORT Cell Fraction")

# grouped boxplot for a chosen cell type (e.g. CD8+ T, Fig 4C)
cib_long <- im_cibersort %>%
  pivot_longer(-c(ID, `P-value_CIBERSORT`, Correlation_CIBERSORT, RMSE_CIBERSORT),
               names_to = "cell_type", values_to = "fraction") %>%
  mutate(cell_type = gsub("_CIBERSORT|_", " ", cell_type)) %>%
  left_join(mrna_sample_filt, by = c("ID" = "sample"))


# ====== 10. GO / KEGG over-representation (up / down) =======================
# helper: run enrichGO + enrichKEGG on a DEG table, split by direction
run_ora <- function(deg, lfc = 0.5) {       # permissive cutoff for enrichment input
  up   <- deg %>% filter(padj < 0.05, log2FoldChange >  lfc) %>% pull(gene) %>% str_to_title()
  down <- deg %>% filter(padj < 0.05, log2FoldChange < -lfc) %>% pull(gene) %>% str_to_title()
  to_entrez <- function(g) bitr(g, "SYMBOL", "ENTREZID", "org.Mm.eg.db")$ENTREZID
  list(
    GO_up   = enrichGO(to_entrez(up),   OrgDb = "org.Mm.eg.db", ont = "all",
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05),
    GO_down = enrichGO(to_entrez(down), OrgDb = "org.Mm.eg.db", ont = "all",
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05),
    KEGG_up   = enrichKEGG(to_entrez(up),   organism = "mmu",
                           pvalueCutoff = 0.05, qvalueCutoff = 1),
    KEGG_down = enrichKEGG(to_entrez(down), organism = "mmu",
                           pvalueCutoff = 0.05, qvalueCutoff = 1)
  )
}
ora_sham_nc   <- run_ora(DEG_sham_nc)
ora_mcao_sham <- run_ora(DEG_mcao_sham)


# ====== 11. GSEA: Hallmark / GO-BP / KEGG (Fig 4D-I) ========================
run_gsea <- function(deg, category, subcategory = NULL,
                     minGS = 10, maxGS = 500) {
  gmt <- msigdbr(species = "Mus musculus",
                 category = category, subcategory = subcategory) %>%
    dplyr::select(gs_name, gene_symbol) %>%
    mutate(gene_symbol = toupper(gene_symbol))
  geneList <- deg$log2FoldChange
  names(geneList) <- toupper(deg$gene)
  geneList <- sort(geneList[!is.na(names(geneList))], decreasing = TRUE)
  GSEA(geneList, TERM2GENE = gmt, minGSSize = minGS, maxGSSize = maxGS,
       pvalueCutoff = 1, pAdjustMethod = "BH", by = "fgsea", verbose = FALSE)
}
gsea_H_mcao_sham    <- run_gsea(DEG_mcao_sham, "H",  NULL, 20, 1000)  # Hallmark
gsea_BP_mcao_sham   <- run_gsea(DEG_mcao_sham, "C5", "BP")            # GO-BP
gsea_KEGG_mcao_sham <- run_gsea(DEG_mcao_sham, "C2", "KEGG")         # KEGG
gsea_H_sham_nc      <- run_gsea(DEG_sham_nc,   "H",  NULL, 20, 1000)

# ridge plot example (Fig 4D/G)
ridgeplot(gsea_H_mcao_sham, showCategory = 20, fill = "pvalue",
          orderBy = "NES", decreasing = FALSE) +
  theme(axis.text.y = element_text(size = 8))
write.xlsx(as.data.frame(gsea_H_mcao_sham@result),
           file.path(dir_res, "GSEA_HALLMARK_mcao_sham.xlsx"))


# ====== 12. Public scRNA-seq GSE264408: CCL5/CCL20 source (Fig 4K, S1 I-J) ==
# Expects a pre-built Seurat object or the 10x RAW + metadata from GEO.
# Minimal processing (matching the manuscript Methods):
process_scRNA <- function(sc) {
  DefaultAssay(sc) <- "RNA"
  sc <- NormalizeData(sc)                                   # LogNormalize
  sc <- FindVariableFeatures(sc, nfeatures = 3000)
  sc <- ScaleData(sc, features = VariableFeatures(sc))
  sc <- RunPCA(sc, features = VariableFeatures(sc))
  sc <- FindNeighbors(sc, dims = 1:30, verbose = FALSE)
  sc <- RunUMAP(sc, dims = 1:30, verbose = FALSE)
  sc
}
# scRNA <- readRDS(file.path(dir_data, "GSE264408_seurat.rds"))
# scRNA <- process_scRNA(scRNA)
# DimPlot(scRNA, group.by = "celltype_subset", reduction = "umap")          # Fig S1 I
# FeaturePlot(scRNA, features = c("Ccl5","Ccl20"), reduction = "umap")      # Fig S1 J
# DotPlot(scRNA, features = c("Ccl5","Ccl20"))                              # Fig 4K


# ====== 13. CIBERSORTx-based ILC deconvolution (Fig S1 C-H) =================
# Signature matrix is built from GSE264408 single-cell reference and the bulk
# rectal TPM matrix is run on the CIBERSORTx web server (S-mode, 100 perms):
#   https://cibersortx.stanford.edu/
# Below: post-processing of the returned fractions.
# TME <- read.csv(file.path(dir_res, "CIBERSORTx_Results.csv"), row.names = 1)
# ilc_cols <- c("ILCreg","ILC1","LTI","ILC3..NCR..","ILC3..NCR...1","ILC2")
# TME$ILC_total <- rowSums(TME[, intersect(ilc_cols, colnames(TME))], na.rm = TRUE)
# TME$sample <- rownames(TME)
# plot_df <- left_join(TME, mrna_sample_filt, by = "sample")
# kruskal.test(ILC_total ~ group, data = plot_df)    # overall test (Fig S1C-H)


# ------------------------------ session info --------------------------------
writeLines(capture.output(sessionInfo()),
           file.path(dir_res, "sessionInfo.txt"))
sessionInfo()
