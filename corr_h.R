library(org.Hs.eg.db)

source("visualization_h.R")

vst_mat  <- readRDS("huntington/data/rds/vst_mat.rds")
metadata <- readRDS("huntington/data/rds/metadata.rds")
res_df   <- readRDS("huntington/data/rds/res_df.rds")

#Computation of correlation
h2afz_entrez <- "3015"
h2afz_expr   <- vst_mat[which(rownames(vst_mat) == h2afz_entrez), ]

ctrl_samples <- colnames(vst_mat)[metadata$condition == "Control"]
hd_samples   <- colnames(vst_mat)[metadata$condition == "HD"]

h2afz_ctrl <- h2afz_expr[ctrl_samples]
h2afz_hd   <- h2afz_expr[hd_samples]

vst_ctrl <- vst_mat[, ctrl_samples]
vst_hd   <- vst_mat[, hd_samples]

# Spearman — Control
cor_ctrl  <- apply(vst_ctrl, 1, function(x)
  cor.test(x, h2afz_ctrl, method = "spearman", exact = FALSE)$estimate)
pval_ctrl <- apply(vst_ctrl, 1, function(x)
  cor.test(x, h2afz_ctrl, method = "spearman", exact = FALSE)$p.value)

# Spearman — HD
cor_hd  <- apply(vst_hd, 1, function(x)
  cor.test(x, h2afz_hd, method = "spearman", exact = FALSE)$estimate)
pval_hd <- apply(vst_hd, 1, function(x)
  cor.test(x, h2afz_hd, method = "spearman", exact = FALSE)$p.value)

# DataFrame
cor_df <- data.frame(
  entrez    = rownames(vst_mat),
  cor_ctrl  = cor_ctrl,
  cor_hd    = cor_hd,
  pval_ctrl = pval_ctrl,
  pval_hd   = pval_hd
)
cor_df$padj_ctrl <- p.adjust(cor_df$pval_ctrl, method = "BH")
cor_df$padj_hd   <- p.adjust(cor_df$pval_hd,   method = "BH")

# Entrez → symbol
symbols       <- mapIds(org.Hs.eg.db,
                        keys     = as.character(cor_df$entrez),
                        column   = "SYMBOL",
                        keytype  = "ENTREZID",
                        multiVals = "first")
cor_df$symbol <- symbols

# grouping
cor_df$group <- "Neither"
cor_df$group[abs(cor_df$cor_ctrl) > 0.4 & cor_df$padj_ctrl < 0.05] <- "Control_only"
cor_df$group[abs(cor_df$cor_hd)   > 0.4 & cor_df$padj_hd   < 0.05] <- "HD_only"
cor_df$group[abs(cor_df$cor_ctrl) > 0.4 & cor_df$padj_ctrl < 0.05 &
               abs(cor_df$cor_hd)   > 0.4 & cor_df$padj_hd   < 0.05] <- "Both"

cat("Control'de anlamlı:", sum(abs(cor_df$cor_ctrl) > 0.4 & cor_df$padj_ctrl < 0.05), "\n")
cat("HD'de anlamlı:",      sum(abs(cor_df$cor_hd)   > 0.4 & cor_df$padj_hd   < 0.05), "\n")
cat("Grup dağılımı:\n"); print(table(cor_df$group))

# save dir
dir.create("huntington/results", recursive = TRUE, showWarnings = FALSE)
write.csv(cor_df, "huntington/results/H2AFZ_correlation_results.csv", row.names = FALSE)

# top 10 genes for HD
hd_specific <- cor_df[cor_df$group == "HD_only", ]
hd_specific <- hd_specific[order(abs(hd_specific$cor_hd), decreasing = TRUE), ]
cat("\nHD'ye özgü top 10 gen:\n")
print(head(hd_specific[, c("symbol", "cor_ctrl", "cor_hd", "padj_hd")], 10))

# gene list for pathway analysis
hd_cor_genes <- cor_df$entrez[cor_df$group %in% c("HD_only", "Both")]
cat("Pathway analizi için gen sayısı:", length(hd_cor_genes), "\n")

# save table
cor_table <- cor_df[cor_df$group %in% c("HD_only", "Both") & !is.na(cor_df$symbol), ]
cor_table <- cor_table[, c("symbol", "entrez", "cor_ctrl", "cor_hd", "padj_ctrl", "padj_hd", "group")]
cor_table <- cor_table[order(abs(cor_table$cor_hd), decreasing = TRUE), ]
write.table(cor_table, "huntington/results/H2AFZ_HD_correlated_genes.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Toplam gen:", nrow(cor_table), "\n")

#save
saveRDS(cor_df,       "huntington/data/rds/cor_df.rds")
saveRDS(hd_cor_genes, "huntington/data/rds/hd_cor_genes.rds")
saveRDS(h2afz_expr,   "huntington/data/rds/h2afz_expr.rds")
