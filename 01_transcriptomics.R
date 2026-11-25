###############################################################################
# 01. DIFFERENTIAL EXPRESSION ANALYSIS - GSE40611
# Parotid gland tissue: 17 pSS, 14 non-pSS sicca, 18 controls
# Affymetrix HG-U133 Plus 2.0 (GPL570), processed in 3 batches
#
#   - RMA normalisation from the 49 raw CEL files (folder data/)
#   - limma with processing batch as a covariate
#   - logFC sign: positive = HIGHER in pSS
#
# Input : data/*.CEL.gz (download GSE40611_RAW.tar from GEO, see README)
# Output: results/DEG/
###############################################################################


# =============================================================================
# 0. PACKAGES AND SETTINGS
# =============================================================================
suppressPackageStartupMessages({
  library(affy)
  library(limma)
  library(GEOquery)
  library(hgu133plus2.db)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(pheatmap)
  library(ggrepel)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Bioconductor packages mask some dplyr verbs; restore them
select <- dplyr::select
first  <- dplyr::first
rename <- dplyr::rename

out_dir <- "results/DEG"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(out_dir, f)

# B-cell genes carried forward to the MR analysis
mr_genes <- c("IGHM", "SELL", "UBD", "MS4A1", "CXCL13", "PAX5", "CR2", "BANK1", "PTPRC")

# DEG thresholds
FDR_CUTOFF   <- 0.05   # Benjamini-Hochberg adjusted P
LOGFC_CUTOFF <- 1      # |log2 fold change| > 1, i.e. > 2-fold change


# =============================================================================
# 1. SAMPLE PHENOTYPES (disease group and processing batch)
# =============================================================================
# The series matrix (~11 MB) is downloaded from GEO once and cached in data/.
# R's default 60-second download timeout is often too short for the NCBI server.
options(timeout = 600)
gse <-getGEO("GSE40611", GSEMatrix = TRUE, getGPL = FALSE, destdir = "data")[[1]]

pheno <- pData(gse) %>%
  transmute(gsm   = geo_accession,
            group = `disease status:ch1`,
            batch = `batch:ch1`)

table(pheno$group)                 # expected: Control 18, pSS 17, Sicca 14
table(pheno$group, pheno$batch)    # groups per batch


# =============================================================================
# 2. READ CEL FILES AND RMA NORMALISATION
# =============================================================================
# RMA = background correction + quantile normalisation + probe-set summarisation (log2)
cel_files <- list.files("data", pattern = "\\.CEL\\.gz$", full.names = TRUE)
stopifnot(length(cel_files) == 49)

raw  <- ReadAffy(filenames = cel_files)
eset <- rma(raw)

# Use GSM IDs as column names and align phenotypes to the expression matrix
colnames(eset) <- sub("_.*", "", colnames(eset))
pheno <- pheno[match(colnames(eset), pheno$gsm), ]
stopifnot(identical(pheno$gsm, colnames(eset)))

group <- factor(pheno$group, levels = c("Control", "Sicca", "pSS"))
batch <- factor(pheno$batch)
expr_all <- exprs(eset)


# =============================================================================
# 3. QUALITY CONTROL
# =============================================================================
group_colours <- c(Control = "#2a78d6", Sicca = "#1baf7a", pSS = "#eb6834")

# 3a. Boxplot: after RMA all arrays should have aligned distributions
png(out("QC_boxplot_RMA.png"), width = 2400, height = 1200, res = 200)
boxplot(expr_all, las = 2, cex.axis = 0.5, outline = FALSE,
        col = group_colours[as.character(group)],
        main = "Expression distribution after RMA normalisation", ylab = "log2 expression")
dev.off()

# 3b. PCA: do samples cluster by disease or by batch?
pca_plot <- function(mat, main_title, file) {
  pc  <- prcomp(t(mat))
  var <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  d   <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], group = group, batch = batch)
  p <- ggplot(d, aes(PC1, PC2, colour = group, shape = batch)) +
    geom_point(size = 3) +
    scale_colour_manual(values = group_colours) +
    labs(x = paste0("PC1 (", var[1], "%)"), y = paste0("PC2 (", var[2], "%)"),
         title = main_title, colour = "Group", shape = "Batch") +
    theme_minimal(base_size = 12)
  ggsave(file, p, width = 7, height = 5.5, dpi = 300, bg = "white")
  p
}
pca_plot(expr_all, "PCA before batch correction", out("QC_PCA_before_batch_correction.png"))

# Batch correction is used ONLY for visualisation (PCA and heatmaps).
# For statistical testing, batch enters the limma model as a covariate (step 5).
expr_vis <- removeBatchEffect(expr_all, batch = batch, design = model.matrix(~ group))
pca_plot(expr_vis, "PCA after batch correction", out("QC_PCA_after_batch_correction.png"))


# =============================================================================
# 4. PROBE FILTERING AND ANNOTATION
# =============================================================================
# 4a. Remove Affymetrix control probes (AFFX-) and lowly expressed probes
#     (keep probes above the lower quartile in at least 14 samples = smallest group)
expr <- expr_all[!grepl("^AFFX", rownames(expr_all)), ]
expr_threshold <- quantile(expr, 0.25)
expr <- expr[rowSums(expr > expr_threshold) >= min(table(group)), ]
c(probes_total = nrow(expr_all), probes_tested = nrow(expr))

# 4b. Probe -> gene symbol and Entrez ID
anno <- AnnotationDbi::select(hgu133plus2.db, keys = rownames(expr),
                              columns = c("SYMBOL", "ENTREZID"), keytype = "PROBEID") %>%
  distinct(PROBEID, .keep_all = TRUE)


# =============================================================================
# 5. DIFFERENTIAL EXPRESSION WITH LIMMA (batch as covariate)
# =============================================================================
design <- model.matrix(~ 0 + group + batch)
colnames(design) <- sub("^group", "", colnames(design))

contrasts <- makeContrasts(
  pSS_vs_Control   = pSS - Control,     # main analysis
  pSS_vs_Sicca     = pSS - Sicca,       # comparator with dryness symptoms but no pSS
  Sicca_vs_Control = Sicca - Control,
  levels = design
)

fit <- lmFit(expr, design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)

# Probe-level table with gene symbols
get_table <- function(coef) {
  topTable(fit, coef = coef, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("PROBEID") %>%
    left_join(anno, by = "PROBEID") %>%
    mutate(status = case_when(
      adj.P.Val < FDR_CUTOFF & logFC >  LOGFC_CUTOFF ~ "Up in pSS",
      adj.P.Val < FDR_CUTOFF & logFC < -LOGFC_CUTOFF ~ "Down in pSS",
      TRUE ~ "Not significant"))
}

# Gene level: one probe per gene (the probe with the smallest P)
collapse_to_gene <- function(tt) {
  tt %>%
    filter(!is.na(SYMBOL)) %>%
    group_by(SYMBOL) %>%
    slice_min(P.Value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(P.Value)
}

deg_probe <- get_table("pSS_vs_Control")
deg_gene  <- collapse_to_gene(deg_probe)

table(deg_gene$status)             # number of up- and down-regulated genes

# MR candidate genes (positive logFC = higher in pSS)
deg_gene %>%
  filter(SYMBOL %in% mr_genes) %>%
  select(SYMBOL, PROBEID, logFC, AveExpr, P.Value, adj.P.Val, status) %>%
  mutate(across(c(logFC, AveExpr), ~ round(.x, 2)),
         across(c(P.Value, adj.P.Val), ~ signif(.x, 3)))

write.csv(deg_probe, out("DEG_pSS_vs_Control_per_probe.csv"), row.names = FALSE)
write.csv(deg_gene,  out("DEG_pSS_vs_Control_per_gene.csv"),  row.names = FALSE)
write.csv(deg_gene %>% filter(status != "Not significant"),
          out("DEG_pSS_vs_Control_significant.csv"), row.names = FALSE)
write.csv(collapse_to_gene(get_table("pSS_vs_Sicca")),
          out("DEG_pSS_vs_Sicca_per_gene.csv"), row.names = FALSE)
write.csv(collapse_to_gene(get_table("Sicca_vs_Control")),
          out("DEG_Sicca_vs_Control_per_gene.csv"), row.names = FALSE)


# =============================================================================
# 6. VOLCANO PLOT
# =============================================================================
status_colours <- c("Up in pSS" = "#eb6834", "Down in pSS" = "#2a78d6",
                    "Not significant" = "grey75")
label_genes <- deg_gene %>%
  filter(status != "Not significant") %>%
  slice_min(P.Value, n = 20) %>%
  bind_rows(deg_gene %>% filter(SYMBOL %in% mr_genes)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

p_volcano <- ggplot(deg_gene, aes(logFC, -log10(P.Value), colour = status)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = "dashed", colour = "grey50") +
  geom_text_repel(data = label_genes, aes(label = SYMBOL), size = 3, colour = "grey15",
                  max.overlaps = 30, min.segment.length = 0) +
  scale_colour_manual(values = status_colours) +
  labs(x = "log2 fold change (pSS vs control)", y = expression(-log[10]~italic(P)), colour = NULL,
       title = "GSE40611: parotid gland, pSS vs control",
       subtitle = paste0("limma adjusted for batch; FDR < ", FDR_CUTOFF,
                         " and |log2FC| > ", LOGFC_CUTOFF)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(out("Volcano_pSS_vs_Control.png"), p_volcano, width = 8, height = 7, dpi = 300, bg = "white")


# =============================================================================
# 7. HEATMAPS (top 20 for the main text, top 50 for the supplement)
# =============================================================================
# Samples ordered by group (Control | Sicca | pSS), with gaps between groups
ord   <- order(group)
annot <- data.frame(Group = group, Batch = batch, row.names = colnames(expr_vis))[ord, ]
gaps  <- cumsum(table(group))[1:2]

make_heatmap <- function(n, main_title, file, width, height, fontsize) {
  top <- deg_gene %>% filter(status != "Not significant") %>% slice_min(P.Value, n = n)
  mat <- expr_vis[top$PROBEID, ord]
  rownames(mat) <- top$SYMBOL
  mat <- t(scale(t(mat)))            # row z-scores
  pheatmap(mat,
           cluster_cols = FALSE, show_colnames = FALSE, gaps_col = gaps,
           annotation_col = annot,
           annotation_colors = list(Group = group_colours,
                                    Batch = c("1" = "grey85", "2" = "grey60", "3" = "grey35")),
           color = colorRampPalette(c("#2a78d6", "white", "#eb6834"))(100),
           breaks = seq(-3, 3, length.out = 101),
           treeheight_row = 25, border_color = NA,
           fontsize = fontsize, fontsize_row = fontsize,
           main = main_title, filename = file, width = width, height = height)
  invisible(top)
}

make_heatmap(20, "Top 20 DEGs (pSS vs control)",
             out("Heatmap_top20.png"), width = 4.6, height = 5.2, fontsize = 9)
make_heatmap(50, "Top 50 DEGs (pSS vs control), row z-scores after batch correction",
             out("Heatmap_top50_supplement.png"), width = 9, height = 10, fontsize = 8)


# =============================================================================
# 8. ENRICHMENT ANALYSIS (GO Biological Process and KEGG)
# =============================================================================
# Background = all tested genes (more appropriate than the whole genome)
universe   <- unique(na.omit(deg_gene$ENTREZID))
genes_up   <- deg_gene %>% filter(status == "Up in pSS")   %>% pull(ENTREZID) %>% na.omit() %>% unique()
genes_down <- deg_gene %>% filter(status == "Down in pSS") %>% pull(ENTREZID) %>% na.omit() %>% unique()
c(up = length(genes_up), down = length(genes_down))

go_up <- enrichGO(genes_up, OrgDb = org.Hs.eg.db, ont = "BP", universe = universe,
                  pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
head(as.data.frame(go_up)[, c("Description", "GeneRatio", "p.adjust")], 15)
write.csv(as.data.frame(go_up), out("GO_BP_up_in_pSS.csv"), row.names = FALSE)

# For the FIGURE only: remove redundant GO terms (simplify) and abbreviate long names.
# The CSV above keeps the full term names.
shorten_go <- function(x) {
  x <- sub("^adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains$",
           "adaptive immune response (Ig-based)", x)
  x <- gsub("positive regulation of", "pos. reg. of", x)
  x <- gsub("negative regulation of", "neg. reg. of", x)
  x <- gsub("regulation of", "reg. of", x)
  x <- gsub("signaling pathway", "signaling", x)
  x <- gsub("immune response-activating", "immune-activating", x)
  x <- gsub("immune response-regulating", "immune-regulating", x)
  x <- gsub("cell surface receptor", "surface receptor", x)
  x <- gsub("antigen processing and presentation", "antigen presentation", x)
  x <- gsub(" involved in immune response", " in immune response", x)
  ifelse(nchar(x) > 45, paste0(substr(x, 1, 42), "..."), x)
}
go_up_plot <- simplify(go_up, cutoff = 0.7, by = "p.adjust", select_fun = min)
go_up_plot@result$Description <- shorten_go(go_up_plot@result$Description)
ggsave(out("GO_BP_up_in_pSS_dotplot.png"),
       dotplot(go_up_plot, showCategory = 15, label_format = 45) +
         ggtitle("GO Biological Process: genes up in pSS"),
       width = 8, height = 7, dpi = 300, bg = "white")

if (length(genes_down) >= 10) {
  go_down <- enrichGO(genes_down, OrgDb = org.Hs.eg.db, ont = "BP", universe = universe,
                      pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
  write.csv(as.data.frame(go_down), out("GO_BP_down_in_pSS.csv"), row.names = FALSE)
}

# KEGG needs internet access to the KEGG server
kegg_up <- tryCatch(
  enrichKEGG(genes_up, organism = "hsa", universe = universe, pAdjustMethod = "BH"),
  error = function(e) { message("KEGG failed: ", conditionMessage(e)); NULL })
if (!is.null(kegg_up)) {
  kegg_up <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(kegg_up), out("KEGG_up_in_pSS.csv"), row.names = FALSE)
  ggsave(out("KEGG_up_in_pSS_dotplot.png"),
         dotplot(kegg_up, showCategory = 15, label_format = 45) +
           ggtitle("KEGG pathways: genes up in pSS"),
         width = 8, height = 7, dpi = 300, bg = "white")
}


# =============================================================================
# 9. SESSION INFORMATION
# =============================================================================
writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
