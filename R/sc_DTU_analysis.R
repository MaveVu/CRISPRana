#' @importFrom SingleCellExperiment counts rowData colData colLabels

# filter cells with exactly 1 gRNA
filter_single_gRNA_cells <- function(sce, gRNA_dir) {
  # Directory contains gRNA files
  files <- list.files(gRNA_dir, pattern = "\\.txt$", full.names = TRUE)

  # Extract gRNAs' names by removing the .txt
  file_labels <- tools::file_path_sans_ext(basename(files))
  cell_lists <- setNames(lapply(files, readLines), file_labels)

  cell_to_files <- list()
  for (grna in names(cell_lists)) {
    cells <- cell_lists[[grna]]
    for (cell in cells) {
      cell_to_files[[cell]] <- c(cell_to_files[[cell]], grna)
    }
  }

  # Identify unique barcodes (appear in exactly one file)
  unique_cells <- names(cell_to_files)[sapply(cell_to_files, length) == 1]

  unique_cell_labels <- sapply(unique_cells, function(cell) cell_to_files[[cell]])

  # Append "-1" to barcodes to match colnames(sce)
  unique_cells_with_suffix <- paste0(unique_cells, "-1")

  labels <- rep(NA_character_, ncol(sce))
  names(labels) <- colnames(sce)

  # Only assign labels to unique barcodes
  for (i in seq_along(unique_cells_with_suffix)) {
    cb <- unique_cells_with_suffix[i]
    if (cb %in% colnames(sce)) {
      labels[cb] <- unique_cell_labels[[i]]
    }
  }

  # Assign to colLabels
  colLabels(sce) <- labels
  # Remove cells with no label
  final_sce <- sce[, !is.na(colLabels(sce))]
  return(final_sce)
}

#' Performing DTU analysis between 2 gRNA groups
#' @param sce The \code{SingleCellExperiment} object generated from the \code{FLAMES} preprocessing step
#' @param gRNA_dir The directory containg gRNA files
#' @param method The chosen statistical method, containing "chisq" (Chi-squared),
#' "multi" (multinomial), and "dir-multi" (Dirichlet-Multinomial)
#' @param gRNA1 The first gRNA of interest
#' @param gRNA2 The second gRNA of interest
#' @param min_count The threshold used in the \code{FLAMES} filtering approach
#' @param filter_by_entry The boolean value to change to \code{CRISPRana} filtering approach
#' @param min_entry The threshold used in the \code{CRISPRana} filtrting approach
#' @return a \code{tibble}
#' @export
sc_DTU_analysis <- function(sce, gRNA_dir, method, gRNA1, gRNA2, min_count = 0.001, filter_by_entry = FALSE, min_entry = 5) {
  # Filter cells with exactly 1 gRNA
  sce <- filter_single_gRNA_cells(sce, gRNA_dir)

  # Filtering step
  if (filter_by_entry == TRUE){
    sub_sce <- sce[, sce$label %in% c(gRNA1, gRNA2)]
    counts_matrix <- counts(sub_sce)
    pseudobulk_matrix <- sapply(unique(colLabels(sub_sce)), function(g) {
      rowSums(counts_matrix[, colLabels(sub_sce) == g, drop = FALSE])
    })
    keep_transcripts <- apply(pseudobulk_matrix, 1, function(row) all(row >= min_entry))
    pseudobulk_matrix <- pseudobulk_matrix[keep_transcripts, ]
    subsce1 <- sub_sce[keep_transcripts, ]
  }else{
    sub_sce <- sce[, sce$label %in% c(gRNA1, gRNA2)]
    counts_matrix <- counts(sub_sce)
    pseudobulk_matrix <- sapply(unique(colLabels(sub_sce)), function(g) {
      rowSums(counts_matrix[, colLabels(sub_sce) == g, drop = FALSE])
    })
    keep_transcripts <- rowSums(pseudobulk_matrix) > min_count
    pseudobulk_matrix <- pseudobulk_matrix[keep_transcripts, ]
    subsce1 <- sub_sce[keep_transcripts, ]
  }

  gene_counts <- table(rowData(subsce1)$gene_id)
  genes_to_keep <- names(gene_counts[gene_counts >= 2])
  keep_rows <- rowData(subsce1)$gene_id %in% genes_to_keep
  subsce1 <- subsce1[keep_rows, ]

  # Chi-squared
  if (method == "chisq"){
    FLAMES::sc_DTU_analysis(subsce1, min_count = 0.001, method = "chisq")
  }
  # Multinomial
  else if (method == "multi"){
    filtered_genes <- unique(rowData(subsce1)$gene_id)
    results <- data.frame(gene_id = character(), statistic = numeric(), df = numeric(), pval = numeric(), stringsAsFactors = FALSE)
    for (gene in filtered_genes) {
      transcripts <- rownames(rowData(subsce1))[rowData(subsce1)$gene_id == gene]
      # round up counts
      counts_mat <- ceiling(counts(subsce1[transcripts, ]))
      df <- as.data.frame(as.matrix(t(counts_mat)))
      df$label <- colLabels(subsce1[transcripts, ])
      df$label <- as.factor(df$label)
      # filter cell with count = 0
      transcript_cols <- setdiff(colnames(df), "label")
      df[transcript_cols] <- lapply(df[transcript_cols], as.numeric)
      df <- df[rowSums(df[transcript_cols]) > 0, ]
      if (nrow(df) < 2) next
      if (length(levels(df$label)) < 2) next

      response_part <- paste0("cbind(", paste(transcript_cols, collapse = ", "), ")")
      formula1 <- as.formula(paste(response_part, "~ label"))
      formula0 <- as.formula(paste(response_part, "~ 1"))

      # Try fitting and testing
      test_res <- tryCatch({
        fit1 <- VGAM::vglm(formula1, family = multinomial(), data = df)
        fit0 <- VGAM::vglm(formula0, family = multinomial(), data = df)
        anova(fit0, fit1, test = "LRT", type = 1)
      }, error = function(e) {
        NULL
      })

      if (!is.null(test_res)) {
        pval <- test_res$`Pr(>Chi)`[2]
        stat <- test_res$`Deviance`[2]
        df <- test_res$Df[2]

        if (!is.na(pval)) {
          results <- rbind(results, data.frame(gene_id = gene, statistic = stat, df = df, pval = pval, stringsAsFactors = FALSE))
        }
      }
    }
    # Adjust p-values using BH
    results$adj_pvalue <- p.adjust(results$pval, method = "BH")
    results <- results[order(results$adj_pvalue), ]
    return(results)
  }
  # Dirichlet-Multinomial
  else if (methods == "dir-multi"){
    counts_df <- as.data.frame(as.matrix(counts(subsce1)), check.names = FALSE)
    counts_df$feature_id <- rowData(subsce1)$transcript_id
    counts_df$gene_id <- rowData(subsce1)$gene_id
    counts_df <- counts_df[, c("gene_id", "feature_id", colnames(counts(subsce1)))]

    samples <- data.frame(
      sample_id = colnames(counts(subsce1)),
      group = colData(subsce1)$label
    )

    d <- DRIMSeq::dmDSdata(counts = counts_df, samples = samples)
    cts <- counts(d)

    design_full <- model.matrix(~ group, data = samples(d_filtered))
    group_name <- colnames(design_full)[2]
    set.seed(123)
    d_filtered <- DRIMSeq::dmPrecision(d_filtered, design = design_full)
    d_filtered <- DRIMSeq::dmFit(d_filtered, design = design_full)
    d_filtered <- DRIMSeq::dmTest(d_filtered, coef=group_name)
    results <- DRIMSeq::results(d_filtered)[order(results(d_filtered)$adj_pvalue), ]
    return(results)
  }
  else{
    stop("Error: Available methods include chisq, multi, dir-multi")
  }
}
