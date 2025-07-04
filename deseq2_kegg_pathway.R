#!/usr/bin/env Rscript

# install.packages("tidyverse", lib = "~/R/x86_64-pc-linux-gnu-library/4.4")
# install.packages(c("tidyverse", "pheatmap", "httr", "jsonlite", "VennDiagram", "UpSetR", "igraph", "ggraph", "gridExtra"), lib = "~/R/x86_64-pc-linux-gnu-library/4.4")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install(version = "3.20")  # Set correct version for R 4.4 BiocManager::install("DESeq2")

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(httr)
  library(jsonlite)
  library(VennDiagram)
  library(UpSetR)
  library(igraph)
  library(ggraph)
  library(gridExtra)
})

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript deseq2_kegg_pathway.R <ko_table.tsv> <metadata.tsv> <output_dir>")
}

ko_file <- args[1]
meta_file <- args[2]
outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# Load data
# -------------------------------
ko <- read.delim(ko_file, row.names = 1, check.names = FALSE)
ko <- ko[!grepl("\\|", rownames(ko)), ]

meta <- read.delim(meta_file)
rownames(meta) <- meta$sample

# Standardize KO column names to match metadata
colnames(ko) <- gsub("_trimmed.*", "", colnames(ko))

# Match sample names
ko <- ko[, intersect(colnames(ko), meta$sample)]
meta <- meta[colnames(ko), ]

ko[] <- lapply(ko, as.numeric)

# -------------------------------
# DESeq2 setup
# -------------------------------
dds <- DESeqDataSetFromMatrix(round(ko), colData = meta, design = ~ condition)
dds <- DESeq(dds)
conditions <- levels(factor(meta$condition))
combinations <- combn(conditions, 2, simplify = FALSE)

# -------------------------------
# KEGG Mapping
# -------------------------------
get_kegg <- function(ko_ids, type = "pathway") {
  out_list <- list()
  for (ko in ko_ids) {
    url <- paste0("http://rest.kegg.jp/link/", type, "/ko:", gsub(":", "", ko))
    r <- try(GET(url), silent = TRUE)
    if (inherits(r, "try-error") || r$status_code != 200) next
    lines <- strsplit(content(r, as = "text"), "\n")[[1]]
    items <- sub(paste0("^", type, ":"), "", sapply(strsplit(lines, "\t"), `[`, 2))
    out_list[[ko]] <- items
  }
  return(out_list)
}

# -------------------------------
# Containers
# -------------------------------
up_kos_list <- list()
down_kos_list <- list()
all_kegg_up <- tibble()
all_kegg_down <- tibble()

# -------------------------------
# Pairwise DE Analysis
# -------------------------------
for (pair in combinations) {
  c1 <- pair[1]; c2 <- pair[2]
  label <- paste0(c2, "_vs_", c1)
  pair_dir <- file.path(outdir, label)
  dir.create(pair_dir, showWarnings = FALSE, recursive = TRUE)

  res <- results(dds, contrast = c("condition", c2, c1))
  res <- res[order(res$padj), ]

  res_up <- res[which(res$padj < 0.05 & res$log2FoldChange > 1), ]
  res_down <- res[which(res$padj < 0.05 & res$log2FoldChange < -1), ]

  write.csv(res, file = file.path(pair_dir, "deseq2_results.csv"))
  write.csv(res_up, file = file.path(pair_dir, "upregulated_kos.csv"))
  write.csv(res_down, file = file.path(pair_dir, "downregulated_kos.csv"))

  up_kos <- rownames(res_up)
  down_kos <- rownames(res_down)
  up_kos_list[[label]] <- up_kos
  down_kos_list[[label]] <- down_kos

  if (length(up_kos) > 0) {
    df_up <- enframe(get_kegg(up_kos, "pathway"), name = "KO", value = "Pathway") %>% unnest(Pathway)
    df_up$Comparison <- label
    all_kegg_up <- bind_rows(all_kegg_up, df_up)
    write.csv(df_up, file = file.path(pair_dir, "ko_to_pathway_up.csv"), row.names = FALSE)

    g_up <- graph_from_data_frame(df_up %>% distinct(KO, Pathway))
    ggsave(file.path(pair_dir, "ko_pathway_tree_up.png"),
           ggraph(g_up, layout = "dendrogram") +
             geom_edge_diagonal() +
             geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
             theme_void(), width = 10, height = 7)

    bubble_up <- df_up %>% count(Pathway, name = "KO_Count")
    ggsave(file.path(pair_dir, "bubble_upregulated_kos.png"),
           ggplot(bubble_up, aes(x = Pathway, y = KO_Count, size = KO_Count, fill = KO_Count)) +
             geom_point(shape = 21, color = "black") +
             theme_minimal() +
             theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
             labs(title = paste("Upregulated KO counts per KEGG Pathway (", label, ")")),
           width = 10, height = 7)

    heat_data_up <- assay(rlog(dds))[up_kos, rownames(meta)[meta$condition %in% pair]]
    heat_data_up <- t(scale(t(heat_data_up)))
    pheatmap(heat_data_up,
             annotation_col = meta[meta$condition %in% pair, "condition", drop = FALSE],
             cluster_rows = TRUE, cluster_cols = TRUE,
             fontsize_row = 7, fontsize_col = 9,
             filename = file.path(pair_dir, "heatmap_upregulated_kos.png"),
             main = paste("Upregulated KOs:", c2, "vs", c1))

    pdf(file.path(pair_dir, "upset_upregulated_kos.pdf"), width = 9, height = 6)
    upset(fromList(list(up_kos)), mainbar.y.label = paste("Upregulated KOs (", label, ")"))
    dev.off()
  }

  if (length(down_kos) > 0) {
    df_down <- enframe(get_kegg(down_kos, "pathway"), name = "KO", value = "Pathway") %>% unnest(Pathway)
    df_down$Comparison <- label
    all_kegg_down <- bind_rows(all_kegg_down, df_down)
    write.csv(df_down, file = file.path(pair_dir, "ko_to_pathway_down.csv"), row.names = FALSE)

    g_down <- graph_from_data_frame(df_down %>% distinct(KO, Pathway))
    ggsave(file.path(pair_dir, "ko_pathway_tree_down.png"),
           ggraph(g_down, layout = "dendrogram") +
             geom_edge_diagonal() +
             geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
             theme_void(), width = 10, height = 7)

    bubble_down <- df_down %>% count(Pathway, name = "KO_Count")
    ggsave(file.path(pair_dir, "bubble_downregulated_kos.png"),
           ggplot(bubble_down, aes(x = Pathway, y = KO_Count, size = KO_Count, fill = KO_Count)) +
             geom_point(shape = 21, color = "black") +
             theme_minimal() +
             theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
             labs(title = paste("Downregulated KO counts per KEGG Pathway (", label, ")")),
           width = 10, height = 7)

    heat_data_down <- assay(rlog(dds))[down_kos, rownames(meta)[meta$condition %in% pair]]
    heat_data_down <- t(scale(t(heat_data_down)))
    pheatmap(heat_data_down,
             annotation_col = meta[meta$condition %in% pair, "condition", drop = FALSE],
             cluster_rows = TRUE, cluster_cols = TRUE,
             fontsize_row = 7, fontsize_col = 9,
             filename = file.path(pair_dir, "heatmap_downregulated_kos.png"),
             main = paste("Downregulated KOs:", c2, "vs", c1))

    pdf(file.path(pair_dir, "upset_downregulated_kos.pdf"), width = 9, height = 6)
    upset(fromList(list(down_kos)), mainbar.y.label = paste("Downregulated KOs (", label, ")"))
    dev.off()
  }
}

# Global Visualizations
visual_dir <- file.path(outdir, "visuals/global")
dir.create(visual_dir, recursive = TRUE, showWarnings = FALSE)

# KO-to-Pathway Trees
if (nrow(all_kegg_up) > 0) {
  g_up_all <- graph_from_data_frame(all_kegg_up %>% distinct(KO, Pathway))
  ggsave(file.path(visual_dir, "global_ko_pathway_tree_up.png"),
         ggraph(g_up_all, layout = "dendrogram") +
           geom_edge_diagonal() +
           geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
           theme_void(), width = 10, height = 7)
}

if (nrow(all_kegg_down) > 0) {
  g_down_all <- graph_from_data_frame(all_kegg_down %>% distinct(KO, Pathway))
  ggsave(file.path(visual_dir, "global_ko_pathway_tree_down.png"),
         ggraph(g_down_all, layout = "dendrogram") +
           geom_edge_diagonal() +
           geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
           theme_void(), width = 10, height = 7)
}

# Bubble Charts
if (nrow(all_kegg_up) > 0) {
  bubble_all_up <- all_kegg_up %>% count(Pathway, Comparison, name = "KO_Count")
  ggsave(file.path(visual_dir, "global_bubble_upregulated_kos.png"),
         ggplot(bubble_all_up, aes(x = Comparison, y = Pathway, size = KO_Count, fill = KO_Count)) +
           geom_point(shape = 21, color = "black") +
           theme_minimal() +
           theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
           labs(title = "Global Upregulated KO counts per KEGG Pathway"),
         width = 10, height = 7)
}

if (nrow(all_kegg_down) > 0) {
  bubble_all_down <- all_kegg_down %>% count(Pathway, Comparison, name = "KO_Count")
  ggsave(file.path(visual_dir, "global_bubble_downregulated_kos.png"),
         ggplot(bubble_all_down, aes(x = Comparison, y = Pathway, size = KO_Count, fill = KO_Count)) +
           geom_point(shape = 21, color = "black") +
           theme_minimal() +
           theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
           labs(title = "Global Downregulated KO counts per KEGG Pathway"),
         width = 10, height = 7)
}

# UpSet Plots
pdf(file.path(visual_dir, "global_upset_upregulated_kos.pdf"), width = 10, height = 7)
upset(fromList(up_kos_list), order.by = "freq", mainbar.y.label = "Global Upregulated KOs")
dev.off()

pdf(file.path(visual_dir, "global_upset_downregulated_kos.pdf"), width = 10, height = 7)
upset(fromList(down_kos_list), order.by = "freq", mainbar.y.label = "Global Downregulated KOs")
dev.off()

# Global Venn Diagrams
venn_dir <- file.path(outdir, "visuals/venn")
dir.create(venn_dir, recursive = TRUE, showWarnings = FALSE)
venn.diagram(x = up_kos_list, filename = file.path(venn_dir, "venn_upregulated.png"), output = TRUE)
venn.diagram(x = down_kos_list, filename = file.path(venn_dir, "venn_downregulated.png"), output = TRUE)

cat("✅ Pairwise and global KO-Pathway visualizations complete.\n")
