#Huntington - H2AZ
#Last Update 27 Feb 2026 - Betül İrem Yardımcı

library(GEOquery)
library(Biobase)
library(org.Hs.eg.db)
library(ggplot2)
library(gridExtra)
library(DESeq2)
library(tidyr)
library(dplyr)

source("visualization_h.R")

histone_entrez_map <- c(H2AFX = "3014", H2AFY = "9555", H2AFY2 = "55506", H2AFZ = "3015")
#histone_entrez_map <- c(H3F3A = "3021", H3F3B = "3022", CENPA = "1058")

OUTPUT_DIR <- "huntington/results/h2a"
PLOT_DIR   <- "huntington/plots/h2a"

#H3 + CENPA için:
#OUTPUT_DIR <- "huntington/results/h3_cenpa"
#PLOT_DIR   <- "huntington/plots/h3_cenpa"


geo_id   <- "GSE64810"
filepath <- paste0("huntington/data/", geo_id, "_raw_counts.tsv")

if (file.exists(filepath)) {
  raw_counts <- read.table(filepath)
} else {
  download.file(
    "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE64810&format=file&file=GSE64810_raw_counts_GRCh38.p13_NCBI.tsv.gz",
    destfile = "GSE64810_raw_counts.tsv.gz",
    mode = "wb"
  )
  raw_counts <- read.table(gzfile("GSE64810_raw_counts.tsv.gz"),
                           header = TRUE, row.names = 1, sep = "\t")
  dir.create("huntington/data", recursive = TRUE, showWarnings = FALSE)
  write.table(raw_counts, filepath, sep = "\t", quote = FALSE)
}

filepath_norm <- paste0("huntington/data/", geo_id, "_norm_counts.tsv")

if (file.exists(filepath_norm)) {
  norm_counts <- read.table(filepath_norm, header = TRUE, sep = "\t")
} else {
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE64nnn/GSE64810/suppl/GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz",
    destfile = "GSE64810_norm_counts.txt.gz",
    mode = "wb"
  )
  norm_counts <- read.table(gzfile("GSE64810_norm_counts.txt.gz"),
                            header = TRUE, row.names = 1, sep = "\t")
  dir.create("huntington/data", recursive = TRUE, showWarnings = FALSE)
  write.table(norm_counts, filepath_norm, sep = "\t", quote = FALSE)
}

# ============================================================
# METADATA
# ============================================================
group    <- ifelse(substr(colnames(norm_counts), 1, 1) == "C", "Control", "HD")
metadata <- data.frame(
  sample    = colnames(raw_counts),
  condition = factor(group, levels = c("Control", "HD"))
)
rownames(metadata) <- colnames(raw_counts)

write.table(metadata,
            paste0("huntington/data/", geo_id, "_metadata.tsv"),
            sep = "\t", quote = FALSE)

rds_path_dds    <- "huntington/data/rds/dds.rds"
rds_path_res    <- "huntington/data/rds/res_df.rds"
rds_path_vst    <- "huntington/data/rds/vst_mat.rds"
rds_path_meta   <- "huntington/data/rds/metadata.rds"

if (file.exists(rds_path_dds)) {
  dds     <- readRDS(rds_path_dds)
  res_df  <- readRDS(rds_path_res)
  vst_mat <- readRDS(rds_path_vst)
  metadata <- readRDS(rds_path_meta)
} else {
  raw_counts_mat <- as.matrix(raw_counts)
  raw_counts_int <- round(raw_counts_mat)
  storage.mode(raw_counts_int) <- "integer"
  
  dds  <- DESeqDataSetFromMatrix(countData = raw_counts_int,
                                 colData   = metadata,
                                 design    = ~ condition)
  keep <- rowSums(counts(dds) >= 10) >= 5
  dds  <- dds[keep, ]
  dds  <- DESeq(dds)
  
  res    <- results(dds, contrast = c("condition", "HD", "Control"), alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$entrez <- rownames(res_df)
  summary(res)
  
  vst_counts <- vst(dds, blind = FALSE)
  vst_mat    <- assay(vst_counts)
  
  dir.create("huntington/data/rds", recursive = TRUE, showWarnings = FALSE)
  saveRDS(dds,      rds_path_dds)
  saveRDS(res_df,   rds_path_res)
  saveRDS(vst_mat,  rds_path_vst)
  saveRDS(metadata, rds_path_meta)
}


histone_results <- res_df[res_df$entrez %in% histone_entrez_map, ]
histone_results$gene <- names(histone_entrez_map)[match(histone_results$entrez, histone_entrez_map)]

cat("\nHistone DESeq2 sonuçları:\n")
print(histone_results[, c("gene", "log2FoldChange", "pvalue", "padj")])

sig_histones <- histone_results[!is.na(histone_results$padj) & histone_results$padj < 0.05, ]
if (nrow(sig_histones) == 0) {
  cat("\nUYARI: Hiçbir histone geni anlamlı değil (padj < 0.05). İleri analiz önerilmez.\n")
} else {
  cat("\nAnlamlı histone genleri:\n")
  print(sig_histones[, c("gene", "log2FoldChange", "padj")])
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(res_df, paste0(OUTPUT_DIR, "/deseq2_full_results.csv"), row.names = FALSE)

#VST EXPRESSION - BAR PLOT
histone_vst <- vst_mat[rownames(vst_mat) %in% histone_entrez_map, , drop = FALSE]
rownames(histone_vst) <- names(histone_entrez_map)[match(rownames(histone_vst), histone_entrez_map)]

histone_long <- as.data.frame(histone_vst) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
  mutate(group = ifelse(substr(sample, 1, 10) %in%
                          colnames(raw_counts)[metadata$condition == "Control"],
                        "Control", "HD"))

histone_stats <- histone_long %>%
  group_by(gene, group) %>%
  summarise(mean_expr = mean(expression),
            sd_expr   = sd(expression),
            n         = n(),
            .groups   = "drop")

sig_data <- res_df %>%
  filter(entrez %in% histone_entrez_map) %>%
  mutate(gene  = names(histone_entrez_map)[match(entrez, histone_entrez_map)]) %>%
  mutate(label = case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE         ~ "ns"
  ))

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

plot_histone_expression(
  histone_stats = histone_stats,
  sig_data      = sig_data,
  gene_list     = names(histone_entrez_map),
  save_path     = paste0(PLOT_DIR, "/histone_expression_barplot.png")
)
