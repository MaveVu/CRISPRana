#' @importFrom SingleCellExperiment counts rowData colData colLabels
generate_sce <- function(final_sce, c1, c2, gene, min_entry = 5){
  # Get the sub-sce for a specific gene and 2 gRNAs
  sub_sce <- final_sce[, colLabels(final_sce) %in% c(c1, c2)]
  gene_id <- subset(gene_map, gene_name == gene)$gene_id
  sub_sce <- sub_sce[rowData(sub_sce)$gene_id == gene_id]
  counts_matrix <- counts(sub_sce)

  # filter trancscripts based on the filtering strategy
  pseudobulk_matrix <- sapply(unique(colLabels(sub_sce)), function(g) {
    rowSums(counts_matrix[, colLabels(sub_sce) == g, drop = FALSE])
  })
  keep_transcripts <- apply(pseudobulk_matrix, 1, function(row) all(row >= min_entry))
  sub_sce <- sub_sce[keep_transcripts, ]
  return (sub_sce)
}

#' Create a heatmap visualizing cell-level transcript usage of a gene between 2 gRNA groups
#' @param final_sce The \code{SingleCellExperiment} object with cell labels
#' @param c1 The first gRNA group
#' @param c2 The second gRNA group
#' @param gene The gene of interest
#' @param min_entry The threshold for the filtering step
#' @param n_clusters The number of clusters for K-means
#' @param random_seed The random seed
#' @return A heatmap showing the cell-level transcript usage of the gene between
#' 2 groups, with a column showing K-means clusters, and a column showing the
#' total molecule counts of the gene in each cell
#' @importFrom circlize colorRamp2 rand_color
#' @importFrom ComplexHeatmap rowAnnotation Heatmap draw
#' @importFrom grid unit grid.grabExpr
#' @importFrom gridExtra grid.arrange
#' @export
plot_heatmap <- function(final_sce, c1, c2, gene, min_entry = 5, n_clusters = 3, random_seed = 0) {
  # Get the sub-sce for a specific gene and 2 gRNAs
  sce <- generate_sce(final_sce, c1, c2, gene, min_entry)

  # Extract transcript count matrix
  expr_mat <- as.matrix(counts(sce))

  # Transpose to get cells x transcripts for normalization
  expr_t <- t(expr_mat)

  # Normalize to proportions
  row_sums <- rowSums(expr_t)
  prop_mat <- sweep(expr_t, 1, row_sums, "/")
  keep_cells <- row_sums > 0
  prop_mat <- prop_mat[keep_cells, ]
  row_sums <- row_sums[keep_cells]

  # Perform k-means clustering on filtered cells using the proportion matrix
  set.seed(random_seed)
  km_result <- kmeans(prop_mat, centers=n_clusters)

  # Store the cluster assignments
  cell_clusters <- km_result$cluster

  # Group labels
  group_labels <- ifelse(colData(sce)$label == c1, c1, c2)[keep_cells]

  # Split matrix by group
  c1_indices <- which(group_labels == c1)
  c2_indices <- which(group_labels == c2)

  c1_mat <- prop_mat[c1_indices, , drop = FALSE]
  c2_mat <- prop_mat[c2_indices, , drop = FALSE]

  # Get clusters for each group
  c1_clusters <- cell_clusters[c1_indices]
  c2_clusters <- cell_clusters[c2_indices]

  c1_counts <- row_sums[c1_indices]
  c2_counts <- row_sums[c2_indices]

  # Sort columns based on overall pattern
  avg_props <- colMeans(prop_mat)
  transcript_order <- order(avg_props, decreasing = TRUE)

  # Apply transcript ordering to both matrices
  c1_mat <- c1_mat[, transcript_order, drop = FALSE]
  c2_mat <- c2_mat[, transcript_order, drop = FALSE]

  # For each cluster, sort cells by their proportion patterns
  c1_sorted_indices <- c()
  for (cluster_id in sort(unique(c1_clusters))) {
    # Get indices of cells in this cluster
    cluster_cells <- which(c1_clusters == cluster_id)

    if (length(cluster_cells) > 1) {
      # Calculate a score for each cell based on weighted transcript usage
      weights <- 1:ncol(c1_mat)
      cell_scores <- rowSums(c1_mat[cluster_cells, ] * matrix(weights,
                                                              nrow=length(cluster_cells),
                                                              ncol=length(weights),
                                                              byrow=TRUE))
      # Order cells by this score
      ordered_indices <- cluster_cells[order(cell_scores)]
      c1_sorted_indices <- c(c1_sorted_indices, ordered_indices)
    } else {
      c1_sorted_indices <- c(c1_sorted_indices, cluster_cells)
    }
  }

  # Repeat for c2
  c2_sorted_indices <- c()
  for (cluster_id in sort(unique(c2_clusters))) {
    cluster_cells <- which(c2_clusters == cluster_id)

    if (length(cluster_cells) > 1) {
      weights <- 1:ncol(c2_mat)
      cell_scores <- rowSums(c2_mat[cluster_cells, ] * matrix(weights,
                                                              nrow=length(cluster_cells),
                                                              ncol=length(weights),
                                                              byrow=TRUE))
      ordered_indices <- cluster_cells[order(cell_scores)]
      c2_sorted_indices <- c(c2_sorted_indices, ordered_indices)
    } else {
      c2_sorted_indices <- c(c2_sorted_indices, cluster_cells)
    }
  }

  # Reorder matrices based on sorted indices
  c1_mat_ordered <- c1_mat[c1_sorted_indices, , drop = FALSE]
  c2_mat_ordered <- c2_mat[c2_sorted_indices, , drop = FALSE]

  c1_counts_ordered <- c1_counts[c1_sorted_indices]
  c2_counts_ordered <- c2_counts[c2_sorted_indices]

  # Get reordered cluster assignments for annotation
  c1_clusters_ordered <- c1_clusters[c1_sorted_indices]
  c2_clusters_ordered <- c2_clusters[c2_sorted_indices]

  # Define color scale from 0 to 1 (since values are proportions)
  col_fun <- colorRamp2(c(0, 1), c("white", "red"))
  # Define color scale for total counts
  count_max <- max(c(c1_counts_ordered, c2_counts_ordered))
  count_min <- min(c(c1_counts_ordered, c2_counts_ordered))
  count_col_fun <- colorRamp2(
    c(count_min, count_max),
    c("lightblue", "darkblue")
  )

  # Create cluster annotation colors
  set.seed(random_seed)
  colors <- rand_color(n_clusters, friendly=TRUE, luminosity="bright")
  # Create named vector of cluster colors
  cluster_ids <- as.character(1:n_clusters)
  cluster_colors <- setNames(colors, cluster_ids)

  # Create annotations
  c1_anno <- rowAnnotation(
    cluster = c1_clusters_ordered,
    col = list(cluster = cluster_colors),
    show_legend = FALSE,
    show_annotation_name = FALSE
  )

  c2_anno <- rowAnnotation(
    cluster = c2_clusters_ordered,
    col = list(cluster = cluster_colors),
    show_legend = FALSE,
    show_annotation_name = FALSE
  )

  c1_count_hm <- Heatmap(
    matrix(c1_counts_ordered, ncol=1),
    name = "Total Count",
    show_row_names = FALSE,
    show_column_names = FALSE,
    width = unit(0.5, "cm"),
    col = count_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  c2_count_hm <- Heatmap(
    matrix(c2_counts_ordered, ncol=1),
    name = "Total Count",
    show_row_names = FALSE,
    show_column_names = FALSE,
    width = unit(0.5, "cm"),
    col = count_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  # Create heatmaps
  ht_c1 <- Heatmap(
    c1_mat_ordered,
    name = "Proportion",
    cluster_rows = FALSE,  # Already clustered and ordered
    cluster_columns = FALSE,  # Already ordered
    show_row_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_title = gene,
    column_title = c1,
    left_annotation = c1_anno,
    col = col_fun,
    use_raster = TRUE
  )

  ht_c2 <- Heatmap(
    c2_mat_ordered,
    name = "Proportion",
    cluster_rows = FALSE,  # Already clustered and ordered
    cluster_columns = FALSE,  # Already ordered
    show_row_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    column_title = c2,
    left_annotation = c2_anno,
    col = col_fun,
    use_raster = TRUE
  )

  ht_list1 <- ht_c1 + c1_count_hm
  ht_list2 <- ht_c2 + c2_count_hm

  # Combine and plot side by side
  g1 <- grid.grabExpr(draw(ht_list1, heatmap_legend_side = "right"))
  g2 <- grid.grabExpr(draw(ht_list2, heatmap_legend_side = "right"))
  grid.arrange(g1, g2, nrow = 1)
}
