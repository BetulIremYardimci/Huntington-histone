#Huntington - H2AZ
#Last Update 27 Feb 2026 - Betül İrem Yardımcı

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(org.Hs.eg.db)
})
histone_genes <- c("H2AX", "H2AZ1", "MACROH2A1", "MACROH2A2")
ensembl_ids <- c(
  "ENSG00000188486",  # H2AX
  "ENSG00000164032",  # H2AZ1
  "ENSG00000134986",  # MACROH2A1
  "ENSG00000172264"   # MACROH2A2
)
gene_map <- setNames(histone_genes, ensembl_ids)

output_dir <- "results/BRCA/variant_coexpression"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load("data/processed/BRCA_se_combat.RData")
se_tumor <- se[, colData(se)$sample_type == "Tumor"]

dds <- DESeqDataSet(se_tumor, design = ~ 1)
dds <- estimateSizeFactors(dds)
vst_obj <- vst(dds, blind = FALSE)
vst_matrix <- assay(vst_obj)

cat("Tumor samples:", ncol(se_tumor), "\n\n")

extract_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) tmp$grobs[[leg]] else NULL
}

analyze_variant_coexpression <- function(target_ensembl, 
                                         target_symbol, 
                                         vst_matrix,
                                         top_n = 10,
                                         cor_threshold = 0.5,
                                         output_dir) {
  
  cat("\n---", target_symbol, "---\n")
  
  # Target expression
  target_expr <- vst_matrix[target_ensembl, ]
  
  # Calculate correlations with ALL genes
  cat("Calculating correlations...\n")
  correlations <- apply(vst_matrix, 1, function(gene_expr) {
    tryCatch({
      cor(target_expr, gene_expr, method = "spearman", use = "complete.obs")
    }, error = function(e) NA)
  })
  
  # Filter and sort
  cor_df <- data.frame(
    ensembl_id = names(correlations),
    correlation = correlations,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(correlation),
           ensembl_id != target_ensembl,
           abs(correlation) > cor_threshold) %>%
    arrange(desc(abs(correlation))) %>%
    head(top_n)
  
  if (nrow(cor_df) == 0) {
    cat("  ⚠ No genes above |r| >", cor_threshold, "\n")
    cat("  Trying lower threshold (0.4)...\n")
    
    cor_df <- data.frame(
      ensembl_id = names(correlations),
      correlation = correlations,
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(correlation),
             ensembl_id != target_ensembl,
             abs(correlation) > 0.4) %>%
      arrange(desc(abs(correlation))) %>%
      head(top_n)
    
    cor_threshold <- 0.4
  }
  
  if (nrow(cor_df) == 0) {
    cat("Still no genes above r > 0.4. Skipping.\n")
    return(NULL)
  }
  
  cat("  Found", nrow(cor_df), "genes above |r| >", cor_threshold, "\n")
  
  # Map to gene symbols
  cor_df$gene_symbol <- mapIds(
    org.Hs.eg.db,
    keys = cor_df$ensembl_id,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # Statistical tests
  cor_df$p_value <- sapply(1:nrow(cor_df), function(i) {
    gene_expr <- vst_matrix[cor_df$ensembl_id[i], ]
    test <- cor.test(target_expr, gene_expr, method = "spearman")
    test$p.value
  })
  
  # Save correlation results
  write.csv(cor_df,
            paste0(output_dir, "/", target_symbol, "_top_correlations.csv"),
            row.names = FALSE)
  
  cat("Saved correlation CSV\n")
  
  # -------------------------
  # SCATTER GRID PLOT
  # -------------------------
  
  cat("  Creating scatter grid...\n")
  
  plot_list <- list()
  
  for (i in 1:nrow(cor_df)) {
    gene_sym <- cor_df$gene_symbol[i]
    gene_ens <- cor_df$ensembl_id[i]
    rho <- cor_df$correlation[i]
    pval <- cor_df$p_value[i]
    
    # Gene expression
    gene_expr <- vst_matrix[gene_ens, ]
    
    # Data frame
    df <- data.frame(
      target_expr = as.numeric(target_expr),
      gene_expr = as.numeric(gene_expr)
    )
    
    # P-value formatting
    p_label <- ifelse(pval < 0.001, "p < 0.001",
                      sprintf("p = %.3f", pval))
    
    # Plot
    p <- ggplot(df, aes(x = target_expr, y = gene_expr)) +
      geom_point(size = 1.5, alpha = 0.4, color = "#2C3E50") +
      geom_smooth(method = "lm", se = TRUE, 
                  color = "#E74C3C", fill = "#E74C3C",
                  linewidth = 1, alpha = 0.2) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      labs(
        title = gene_sym,
        subtitle = sprintf("r = %.3f | %s", rho, p_label),
        x = paste0(target_symbol, " (VST)"),
        y = paste0(gene_sym, " (VST)")
      )
    
    plot_list[[gene_sym]] <- p
  }
  n <- nrow(cor_df)
  ncols <- ifelse(n <= 3, n, 3)
  nrows <- ceiling(n / ncols)
  
  # Title
  title_grob <- textGrob(
    paste0(target_symbol, " Co-expression Network in BRCA"),
    gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  subtitle_grob <- textGrob(
    paste0("Top ", n, " correlated genes (|r| > ", cor_threshold, 
           ", n=1090 tumors)"),
    gp = gpar(fontsize = 11, col = "gray40")
  )
  
  # Grid
  plot_grid <- do.call(arrangeGrob, c(plot_list, list(ncol = ncols)))
  
  plot_height <- 4 * nrows + 2
  
  final_plot <- arrangeGrob(
    title_grob,
    subtitle_grob,
    plot_grid,
    heights = c(0.5, 0.3, nrows * 4)
  )
  
  # Save
  save_path <- paste0(output_dir, "/", target_symbol, "_scatter_grid.pdf")
  ggsave(save_path, final_plot, 
         width = 12, height = plot_height, 
         dpi = 300, bg = "white")
  
  cat("  ✓ Saved scatter grid:", save_path, "\n")
  
  return(cor_df)
}

all_results <- list()

for (i in 1:length(ensembl_ids)) {
  ens_id <- ensembl_ids[i]
  gene_name <- gene_map[ens_id]
  
  result <- analyze_variant_coexpression(
    target_ensembl = ens_id,
    target_symbol = gene_name,
    vst_matrix = vst_matrix,
    top_n = 10,
    cor_threshold = 0.5,
    output_dir = output_dir
  )
  
  if (!is.null(result)) {
    all_results[[gene_name]] <- result
  }
}

summary_df <- data.frame(
  Variant = character(),
  N_Correlated = integer(),
  Top_Gene = character(),
  Top_Correlation = numeric(),
  stringsAsFactors = FALSE
)

for (gene_name in names(all_results)) {
  res <- all_results[[gene_name]]
  
  summary_df <- rbind(summary_df, data.frame(
    Variant = gene_name,
    N_Correlated = nrow(res),
    Top_Gene = res$gene_symbol[1],
    Top_Correlation = round(res$correlation[1], 3)
  ))
}

print(summary_df)

write.csv(summary_df,
          paste0(output_dir, "/summary_all_variants.csv"),
          row.names = FALSE)


library(clusterProfiler)
library(org.Hs.eg.db)

h2ax_top_genes <- c("UBE2S", "SAC3D1", "CDT1", "SNRPA", "CREBRF", 
                    "HMBS", "POP7", "PKMYT1", "UBE2C", "PPP1R14B")

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = h2ax_top_genes,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

kegg_result <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

go_result <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

print(kegg_result@result[, c("Description", "p.adjust", "geneID")])

print(go_result@result[, c("Description", "p.adjust", "geneID")])


library(enrichplot)

dotplot(go_result, showCategory = 15) +
  labs(title = "H2AX Co-expression Network",
       subtitle = "GO Biological Process Enrichment")

ggsave("results/BRCA/h2ax_network/pathway_enrichment_GO.pdf",
       width = 10, height = 8)

dotplot(kegg_result, showCategory = 10) +
  labs(title = "H2AX Co-expression Network",
       subtitle = "KEGG Pathway Enrichment")

ggsave("results/BRCA/h2ax_network/pathway_enrichment_KEGG.pdf",
       width = 10, height = 6)
