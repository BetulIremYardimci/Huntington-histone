library(GEOquery)
library(Biobase)
library(org.Hs.eg.db)
library(ggplot2)
library(gridExtra)
library(DESeq2)
library(tidyr)
library(dplyr)

source("visualization_h.R")

geo_id    <- "GSE64810"
filepath  <- paste0("huntington/data/", geo_id, "_raw_counts.tsv")

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

# METADATA
group    <- ifelse(substr(colnames(norm_counts), 1, 1) == "C", "Control", "HD")
metadata <- data.frame(
  sample    = colnames(raw_counts),
  condition = factor(group, levels = c("Control", "HD"))
)
rownames(metadata) <- colnames(raw_counts)

write.table(metadata,
            paste0("huntington/data/", geo_id, "_metadata.tsv"),
            sep = "\t", quote = FALSE)

# DESEQ2
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

# Histone results
histone_entrez_map <- c(H2AFX = "3014", H2AFY = "9555", H2AFY2 = "55506", H2AFZ = "3015")

histone_results <- res_df[res_df$entrez %in% histone_entrez_map, ]
histone_results$gene <- names(histone_entrez_map)[match(histone_results$entrez, histone_entrez_map)]
print(histone_results[, c("gene", "log2FoldChange", "pvalue", "padj")])

write.csv(res_df, "huntington/results/deseq2_full_results.csv")

# VST
vst_counts <- vst(dds, blind = FALSE)
vst_mat    <- assay(vst_counts)


# Histone VST Expression
histone_vst <- vst_mat[rownames(vst_mat) %in% histone_entrez_map, ]
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

#plot
plot_histone_expression(
  histone_stats = histone_stats,
  sig_data      = sig_data,
  save_path     = "huntington/plots/histone_expression_barplot.png"
)


#save
dir.create("huntington/data/rds", recursive = TRUE, showWarnings = FALSE)
saveRDS(vst_mat,  "huntington/data/rds/vst_mat.rds")
saveRDS(metadata, "huntington/data/rds/metadata.rds")
saveRDS(res_df,   "huntington/data/rds/res_df.rds")
saveRDS(dds,      "huntington/data/rds/dds.rds")
