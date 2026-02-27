library(ggplot2)
library(gridExtra)

COLORS <- c("Control" = "#2980B9", "HD" = "#C0392B")


#1- BAR PLOT
plot_histone_expression <- function(histone_stats, sig_data,
                                    save_path = "huntington/plots/histone_expression_barplot.png") {

  gene_list  <- c("H2AFX", "H2AFY", "H2AFY2", "H2AFZ")
  plot_list  <- list()

  for (g in gene_list) {
    df_plot   <- histone_stats[histone_stats$gene == g, ]
    sig_label <- sig_data$label[sig_data$gene == g]
    p_val_raw <- sig_data$padj[sig_data$gene == g]

    y_max     <- max(df_plot$mean_expr + df_plot$sd_expr)
    y_bracket <- y_max * 1.05
    y_star    <- y_max * 1.12

    p <- ggplot(df_plot, aes(x = group, y = mean_expr, fill = group)) +
      geom_bar(stat = "identity", width = 0.55, alpha = 0.9,
               color = "white", linewidth = 0.5) +
      geom_errorbar(aes(ymin = mean_expr - sd_expr,
                        ymax = mean_expr + sd_expr),
                    width = 0.15, linewidth = 0.7, color = "gray30") +
      annotate("segment", x = 1, xend = 2,
               y = y_bracket, yend = y_bracket, linewidth = 0.6) +
      annotate("segment", x = 1, xend = 1,
               y = y_bracket * 0.98, yend = y_bracket, linewidth = 0.6) +
      annotate("segment", x = 2, xend = 2,
               y = y_bracket * 0.98, yend = y_bracket, linewidth = 0.6) +
      annotate("text", x = 1.5, y = y_star,
               label = sig_label, size = 6, fontface = "bold") +
      scale_fill_manual(values = c("Control" = "#A8C4D4", "HD" = "#F2A8A8")) +
      labs(subtitle = paste0("adj.p = ", format.pval(p_val_raw, digits = 2)),
           y = "VST Normalized Expression") +
      theme_classic() +
      theme(
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30"),
        axis.title.x  = element_blank(),
        legend.position = "none"
      ) +
      ggtitle(g)

    plot_list[[g]] <- p
  }

  title_grob    <- grid::textGrob("Histone Variant Expression — HD vs Control",
                                  gp = grid::gpar(fontsize = 14, fontface = "bold"))
  subtitle_grob <- grid::textGrob("VST normalized counts, DESeq2 adjusted p-values",
                                  gp = grid::gpar(fontsize = 10, col = "gray40"))

  legend_plot <- ggplot(data.frame(group = c("Control", "HD"), x = 1, y = 1),
                        aes(x = x, y = y, fill = group)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Control" = "#A8C4D4", "HD" = "#F2A8A8"), name = "") +
    theme_void() +
    theme(legend.position = "bottom", legend.text = element_text(size = 11))
  legend_grob <- .extract_legend(legend_plot)

  final_plot <- gridExtra::arrangeGrob(
    title_grob, subtitle_grob,
    gridExtra::arrangeGrob(plot_list[["H2AFX"]], plot_list[["H2AFY"]], ncol = 2),
    gridExtra::arrangeGrob(plot_list[["H2AFY2"]], plot_list[["H2AFZ"]], ncol = 2),
    legend_grob,
    heights = c(0.06, 0.04, 1, 1, 0.12)
  )

  dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(save_path, final_plot, width = 12, height = 10, dpi = 300)
  cat("Kaydedildi:", save_path, "\n")
}


#2- KEGG DOT PLOT
plot_kegg_dotplot <- function(kegg_results, top_n = 15,
                              save_path = "huntington/plots/kegg_dotplot.png") {

  kegg_df <- kegg_results@result[kegg_results@result$p.adjust < 0.05, ]
  kegg_df$GeneRatio_num <- sapply(kegg_df$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  })

  kegg_top <- head(kegg_df[order(kegg_df$p.adjust), ], top_n)
  kegg_top$Description <- factor(kegg_top$Description, levels = rev(kegg_top$Description))

  p <- ggplot(kegg_top, aes(x = GeneRatio_num, y = Description,
                            size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#C0392B", high = "#2980B9", name = "p.adjust") +
    scale_size_continuous(range = c(3, 10), name = "Gene Count") +
    theme_classic() +
    theme(axis.text.y  = element_text(size = 10),
          plot.title   = element_text(size = 13, face = "bold")) +
    labs(title = "KEGG Pathway — H2AFZ Correlated Genes (HD)",
         x = "Gene Ratio", y = NULL)

  dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(save_path, p, width = 10, height = 7, dpi = 300)
  cat("Kaydedildi:", save_path, "\n")

  invisible(p)
}


#3- GO DOT PLOT
plot_go_dotplot <- function(go_results, top_n = 15,
                            save_path = "huntington/plots/go_dotplot.png") {

  go_df <- go_results@result[go_results@result$p.adjust < 0.05, ]
  go_df$GeneRatio_num <- sapply(go_df$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  })

  go_top <- head(go_df[order(go_df$p.adjust), ], top_n)
  go_top$Description <- factor(go_top$Description, levels = rev(go_top$Description))

  p <- ggplot(go_top, aes(x = GeneRatio_num, y = Description,
                          size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#C0392B", high = "#2980B9", name = "p.adjust") +
    scale_size_continuous(range = c(3, 10), name = "Gene Count") +
    theme_classic() +
    theme(axis.text.y  = element_text(size = 10),
          plot.title   = element_text(size = 13, face = "bold")) +
    labs(title = "GO Biological Process — H2AFZ Correlated Genes (HD)",
         x = "Gene Ratio", y = NULL)

  dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(save_path, p, width = 11, height = 7, dpi = 300)
  cat("Kaydedildi:", save_path, "\n")

  invisible(p)
}


#4- MITO SCATTER GRID
plot_mito_scatter_grid <- function(target_genes, vst_mat, h2afz_expr, metadata,
                                   save_path = "huntington/plots/H2AFZ_mito_scatter_grid.png") {

  plot_list <- list()

  for (i in 1:nrow(target_genes)) {
    gene_sym    <- target_genes$symbol[i]
    gene_entrez <- as.character(target_genes$entrez[i])
    rho_ctrl    <- round(target_genes$cor_ctrl[i], 2)
    rho_hd      <- round(target_genes$cor_hd[i], 2)

    gene_idx  <- which(rownames(vst_mat) == gene_entrez)
    gene_expr <- as.numeric(vst_mat[gene_idx, ])

    df <- data.frame(
      H2AFZ     = as.numeric(h2afz_expr),
      gene_expr = gene_expr,
      group     = metadata$condition
    )

    p <- ggplot(df, aes(x = H2AFZ, y = gene_expr, color = group)) +
      geom_point(size = 1.8, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
      scale_color_manual(values = COLORS) +
      theme_classic() +
      theme(
        plot.title      = element_text(size = 10, face = "bold"),
        plot.subtitle   = element_text(size = 8, color = "gray40"),
        axis.title      = element_text(size = 8),
        axis.text       = element_text(size = 7),
        legend.position = "none"
      ) +
      labs(
        title    = gene_sym,
        subtitle = paste0("Control: ρ=", rho_ctrl, " | HD: ρ=", rho_hd),
        x        = "H2AFZ (VST)",
        y        = paste0(gene_sym, " (VST)")
      )

    plot_list[[gene_sym]] <- p
  }

  # Legend
  legend_plot <- ggplot(data.frame(group = c("Control", "HD"), x = 1, y = 1),
                        aes(x = x, y = y, color = group)) +
    geom_point() +
    scale_color_manual(values = COLORS, name = "") +
    theme_void() +
    theme(legend.position = "bottom", legend.text = element_text(size = 11))
  legend_grob <- .extract_legend(legend_plot)

  title_grob <- grid::textGrob(
    "H2AFZ vs Mitochondrial Genes — HD-Specific Correlations",
    gp = grid::gpar(fontsize = 13, fontface = "bold")
  )
  subtitle_grob <- grid::textGrob(
    "Genes correlated with H2AFZ exclusively in HD (GO: mitochondrial respirasome assembly)",
    gp = grid::gpar(fontsize = 10, col = "gray40")
  )

  n <- nrow(target_genes)
  grob_list <- plot_list[1:min(n, 9)]
  # 9'dan az gen varsa boş plot ile doldur
  while (length(grob_list) < 9) {
    grob_list[[length(grob_list) + 1]] <- ggplot() + theme_void()
  }

  final_plot <- gridExtra::arrangeGrob(
    title_grob,
    subtitle_grob,
    gridExtra::arrangeGrob(
      grob_list[[1]], grob_list[[2]], grob_list[[3]],
      grob_list[[4]], grob_list[[5]], grob_list[[6]],
      grob_list[[7]], grob_list[[8]], grob_list[[9]],
      ncol = 3
    ),
    legend_grob,
    heights = c(0.08, 0.05, 1, 0.08)
  )

  dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(save_path, final_plot, width = 12, height = 13, dpi = 300)
  cat("Kaydedildi:", save_path, "\n")
}


#helper function
.extract_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}
