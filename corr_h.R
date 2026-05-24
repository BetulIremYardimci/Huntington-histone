#Huntington - H2AZ
#Last Update 27 Feb 2026 - Betül İrem Yardımcı

library(org.Hs.eg.db)
library(AnnotationDbi)

source("visualization_h.R")

TARGET_SYMBOL <- "H2AFZ" 
TARGET_ENTREZ <- "3015"   

# TARGET_SYMBOL <- "H3F3A" ; TARGET_ENTREZ <- "3021"
# TARGET_SYMBOL <- "H3F3B" ; TARGET_ENTREZ <- "3022"
# TARGET_SYMBOL <- "CENPA" ; TARGET_ENTREZ <- "1058"

vst_mat  <- readRDS("huntington/data/rds/vst_mat.rds")
metadata <- readRDS("huntington/data/rds/metadata.rds")
res_df   <- readRDS("huntington/data/rds/res_df.rds")

target_expr <- vst_mat[which(rownames(vst_mat) == TARGET_ENTREZ), ]

ctrl_samples <- colnames(vst_mat)[metadata$condition == "Control"]
hd_samples   <- colnames(vst_mat)[metadata$condition == "HD"]

target_ctrl <- target_expr[ctrl_samples]
target_hd   <- target_expr[hd_samples]

vst_ctrl <- vst_mat[, ctrl_samples]
vst_hd   <- vst_mat[, hd_samples]

# Spearman — Control
cor_ctrl  <- apply(vst_ctrl, 1, function(x)
  cor.test(x, target_ctrl, method = "spearman", exact = FALSE)$estimate)
pval_ctrl <- apply(vst_ctrl, 1, function(x)
  cor.test(x, target_ctrl, method = "spearman", exact = FALSE)$p.value)

# Spearman — HD
cor_hd  <- apply(vst_hd, 1, function(x)
  cor.test(x, target_hd, method = "spearman", exact = FALSE)$estimate)
pval_hd <- apply(vst_hd, 1, function(x)
  cor.test(x, target_hd, method = "spearman", exact = FALSE)$p.value)

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
                        keys      = as.character(cor_df$entrez),
                        column    = "SYMBOL",
                        keytype   = "ENTREZID",
                        multiVals = "first")
cor_df$symbol <- symbols

# Grup ataması
cor_df$group <- "Neither"
cor_df$group[abs(cor_df$cor_ctrl) > 0.4 & cor_df$padj_ctrl < 0.05] <- "Control_only"
cor_df$group[abs(cor_df$cor_hd)   > 0.4 & cor_df$padj_hd   < 0.05] <- "HD_only"
cor_df$group[abs(cor_df$cor_ctrl) > 0.4 & cor_df$padj_ctrl < 0.05 &
               abs(cor_df$cor_hd) > 0.4 & cor_df$padj_hd   < 0.05] <- "Both"

cat("Control'de anlamlı:", sum(abs(cor_df$cor_ctrl) > 0.4 & cor_df$padj_ctrl < 0.05), "\n")
cat("HD'de anlamlı:",      sum(abs(cor_df$cor_hd)   > 0.4 & cor_df$padj_hd   < 0.05), "\n")
cat("Grup dağılımı:\n");   print(table(cor_df$group))

dir.create("huntington/results", recursive = TRUE, showWarnings = FALSE)

write.csv(cor_df,
          paste0("huntington/results/", TARGET_SYMBOL, "_correlation_results.csv"),
          row.names = FALSE)

# HD top 10
hd_specific <- cor_df[cor_df$group == "HD_only", ]
hd_specific <- hd_specific[order(abs(hd_specific$cor_hd), decreasing = TRUE), ]
cat("\nHD'ye özgü top 10 gen:\n")
print(head(hd_specific[, c("symbol", "cor_ctrl", "cor_hd", "padj_hd")], 10))

# Pathway için gen listesi
hd_cor_genes <- cor_df$entrez[cor_df$group %in% c("HD_only", "Both")]
cat("Pathway analizi için gen sayısı:", length(hd_cor_genes), "\n")

# Tablo
cor_table <- cor_df[cor_df$group %in% c("HD_only", "Both") & !is.na(cor_df$symbol), ]
cor_table <- cor_table[, c("symbol", "entrez", "cor_ctrl", "cor_hd", "padj_ctrl", "padj_hd", "group")]
cor_table <- cor_table[order(abs(cor_table$cor_hd), decreasing = TRUE), ]
write.table(cor_table,
            paste0("huntington/results/", TARGET_SYMBOL, "_HD_correlated_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Toplam gen:", nrow(cor_table), "\n")


dir.create("huntington/data/rds", recursive = TRUE, showWarnings = FALSE)
saveRDS(cor_df,       paste0("huntington/data/rds/", TARGET_SYMBOL, "_cor_df.rds"))
saveRDS(hd_cor_genes, paste0("huntington/data/rds/", TARGET_SYMBOL, "_hd_cor_genes.rds"))
saveRDS(target_expr,  paste0("huntington/data/rds/", TARGET_SYMBOL, "_expr.rds"))
cat("RDS dosyaları kaydedildi.\n")