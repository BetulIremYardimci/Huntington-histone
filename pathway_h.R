library(clusterProfiler)
library(org.Hs.eg.db)

source("visualization_h.R")

cor_df       <- readRDS("huntington/data/rds/cor_df.rds")
hd_cor_genes <- readRDS("huntington/data/rds/hd_cor_genes.rds")
vst_mat      <- readRDS("huntington/data/rds/vst_mat.rds")
metadata     <- readRDS("huntington/data/rds/metadata.rds")
h2afz_expr   <- readRDS("huntington/data/rds/h2afz_expr.rds")

# PATHWAY Analysis
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

cat("GO significant pathway:",   nrow(go_results@result[go_results@result$p.adjust < 0.05, ]), "\n")
cat("KEGG significant pathway:", nrow(kegg_results@result[kegg_results@result$p.adjust < 0.05, ]), "\n")
cat("\nTop 10 GO:\n");   print(head(go_results@result[, c("Description", "GeneRatio", "p.adjust")], 10))
cat("\nTop 10 KEGG:\n"); print(head(kegg_results@result[, c("Description", "GeneRatio", "p.adjust")], 10))

#PLOT - PATHWAY
plot_kegg_dotplot(kegg_results, save_path = "huntington/plots/kegg_dotplot.png")
plot_go_dotplot(go_results,     save_path = "huntington/plots/go_dotplot.png")

# PATHWAY × CORRELATION
hd_pathway_genes <- kegg_results@result["hsa05016", "geneID"]
hd_genes_list    <- strsplit(hd_pathway_genes, "/")[[1]]

mito_pathway_genes <- go_results@result["GO:0097250", "geneID"]
mito_genes_list    <- strsplit(mito_pathway_genes, "/")[[1]]

hd_specific_genes <- cor_df[cor_df$group == "HD_only" & !is.na(cor_df$symbol), ]
hd_specific_genes <- hd_specific_genes[order(abs(hd_specific_genes$cor_hd), decreasing = TRUE), ]

overlap_hd <- hd_specific_genes[hd_specific_genes$symbol %in% hd_genes_list, ]
cat("HD-specific + Huntington pathway:", nrow(overlap_hd), "\n")
print(head(overlap_hd[, c("symbol", "cor_ctrl", "cor_hd")], 10))

overlap_mito <- hd_specific_genes[hd_specific_genes$symbol %in% mito_genes_list, ]
cat("HD-specific + Mito pathway:", nrow(overlap_mito), "\n")
print(head(overlap_mito[, c("symbol", "cor_ctrl", "cor_hd")], 10))

#plot
# MİTO SCATTER GRID
target_genes <- overlap_mito[, c("symbol", "entrez", "cor_ctrl", "cor_hd")]

plot_mito_scatter_grid(
  target_genes = target_genes,
  vst_mat      = vst_mat,
  h2afz_expr   = h2afz_expr,
  metadata     = metadata,
  save_path    = "huntington/plots/H2AFZ_mito_scatter_grid.png"
)
