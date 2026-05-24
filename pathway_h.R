#Huntington - H2AZ
#Last Update 27 Feb 2026 - Betül İrem Yardımcı

library(clusterProfiler)
library(org.Hs.eg.db)

source("visualization_h.R")

TARGET_SYMBOL <- "H2AFZ"
TARGET_ENTREZ <- "3015"

# TARGET_SYMBOL <- "H3F3A" ; TARGET_ENTREZ <- "3021"
# TARGET_SYMBOL <- "H3F3B" ; TARGET_ENTREZ <- "3022"
# TARGET_SYMBOL <- "CENPA" ; TARGET_ENTREZ <- "1058"

cor_df       <- readRDS(paste0("huntington/data/rds/", TARGET_SYMBOL, "_cor_df.rds"))
hd_cor_genes <- readRDS(paste0("huntington/data/rds/", TARGET_SYMBOL, "_hd_cor_genes.rds"))
target_expr  <- readRDS(paste0("huntington/data/rds/", TARGET_SYMBOL, "_expr.rds"))
vst_mat      <- readRDS("huntington/data/rds/vst_mat.rds")
metadata     <- readRDS("huntington/data/rds/metadata.rds")


go_results <- enrichGO(
  gene          = as.character(hd_cor_genes),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

kegg_results <- enrichKEGG(
  gene          = as.character(hd_cor_genes),
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

cat("GO anlamlı pathway:",   nrow(go_results@result[go_results@result$p.adjust < 0.05, ]), "\n")
cat("KEGG anlamlı pathway:", nrow(kegg_results@result[kegg_results@result$p.adjust < 0.05, ]), "\n")
cat("\nTop 10 GO:\n");   print(head(go_results@result[, c("Description", "GeneRatio", "p.adjust")], 10))
cat("\nTop 10 KEGG:\n"); print(head(kegg_results@result[, c("Description", "GeneRatio", "p.adjust")], 10))


plot_kegg_dotplot(kegg_results,
                  save_path = paste0("huntington/plots/", TARGET_SYMBOL, "_kegg_dotplot.png"))
plot_go_dotplot(go_results,
                save_path = paste0("huntington/plots/", TARGET_SYMBOL, "_go_dotplot.png"))


# TOP 10 GO PATHWAY × CORRELATION-  SCATTER GRID
go_sig <- go_results@result[go_results@result$p.adjust < 0.05, ]
go_sig$rank_score <- rank(go_sig$p.adjust) - rank(go_sig$Count)
go_top10 <- head(go_sig[order(go_sig$rank_score), ], 10)

cat("Seçilen top 10 GO pathway:\n")
print(go_top10[, c("Description", "Count", "p.adjust")])

dir.create(paste0("huntington/plots/", TARGET_SYMBOL, "_go_scatter_grids"),
           recursive = TRUE, showWarnings = FALSE)

for (i in 1:nrow(go_top10)) {
  
  go_id   <- rownames(go_top10)[i]
  go_name <- go_top10$Description[i]
  
  pathway_genes_str <- go_top10[i, "geneID"]
  if (is.na(pathway_genes_str) || pathway_genes_str == "") next
  pathway_gene_symbols <- strsplit(pathway_genes_str, "/")[[1]]
  
  hd_specific_genes <- cor_df[cor_df$group == "HD_only" & !is.na(cor_df$symbol), ]
  overlap <- hd_specific_genes[hd_specific_genes$symbol %in% pathway_gene_symbols, ]
  
  if (nrow(overlap) == 0) {
    cat("Örtüşen gen yok, atlanıyor:", go_name, "\n")
    next
  }
  
  overlap      <- overlap[order(abs(overlap$cor_hd), decreasing = TRUE), ]
  target_genes <- overlap[1:min(9, nrow(overlap)), c("symbol", "entrez", "cor_ctrl", "cor_hd")]
  
  safe_name <- gsub("[^a-zA-Z0-9_]", "_", go_name)
  save_path <- paste0("huntington/plots/", TARGET_SYMBOL, "_go_scatter_grids/",
                      go_id, "_", safe_name, ".png")
  
  cat("Çiziliyor:", go_name, "(", nrow(overlap), "gen )\n")
  
  plot_mito_scatter_grid(
    target_genes = target_genes,
    vst_mat      = vst_mat,
    h2afz_expr   = target_expr,
    metadata     = metadata,
    pathway_name = go_name,
    save_path    = save_path
  )
}

