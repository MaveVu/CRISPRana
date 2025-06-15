#' @importFrom magrittr %>%
rev_comp <- function(seq) {
  chartr("ATGC", "TACG", seq) %>%
    strsplit("") %>%
    unlist() %>%
    rev() %>%
    paste(collapse = "")
}

#' Generaate gRNA files containing gRNA-cell barcode mapping information
#' @param ref_sgRNA The reference file containing gRNA names and sequences
#' @param gRNA_R1 The Read 1 data containing cell barcodes and UMIs
#' @param gRNA_R2 The Read 2 data containing gRNA sequences
#' @param CB_file The cell barcode whitelist file
#' @param cleared_CB The file containing cleared cell barcodes (no "-1")
#' @param num_counts The threshold for filtering low-counts cell barcodes
#' @param outdir The directory containing gRNA files
#' @return A directory containing gRNA files
#' @export
gRNA_assign <- function(ref_sgRNA, gRNA_R1, gRNA_R2, CB_file, cleared_CB, num_counts = 10, outdir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }

  # Read sgRNA reference file
  sgRNA_data <- read.table(ref_sgRNA, header=FALSE, sep="\t")
  gene_list <- sgRNA_data$V1
  sg_list <- sgRNA_data$V2


  # Step 1
  cat("Step 1: Capturing cell barcodes + UMIs from paired-end gRNA sequencing data...\n")
  R1 <- readLines(gRNA_R1)
  R2 <- readLines(gRNA_R2)
  # Process each gRNA
  all_barcodes <- list()
  for (i in 1:length(sg_list)) {
    sgRNA <- sg_list[i]
    sgRNA_rc <- rev_comp(sgRNA)
    barcodes <- c()

    # FASTQ format
    for (j in 0:(length(R1)/4 - 1)) {
      R1_seq <- R1[4*j + 2]
      R2_seq <- R2[4*j + 2]

      # Check if gRNA reverse complement in R2
      if (grepl(sgRNA_rc, R2_seq, fixed = TRUE)) {
        # Extract first 26 bases from R1
        barcode <- substr(R1_seq, 1, 26)
        barcodes <- c(barcodes, barcode)
      }
    }
    # store barcodes for this gRNA
    all_barcodes[[sgRNA]] <- barcodes
  }

  # Step 2
  cat("Step 2: Consolidating UMIs...\n")
  unique_sg_barcodes <- lapply(all_barcodes, unique)

  # Step 3
  cat("Step 3: Filtering barcodes based on the cell barcode whitelist...\n")
  # Check if cleared barcode file exists. If not, create it from the original cell barcode whitelist
  if (!file.exists(cleared_CB)) {
    cat("Cleared barcode file doesn't exist. Creating it from the original cell barcode whitelist...\n")
    if (!file.exists(CB_file)) {
      stop("Error: The cell barcode whitelist file doesn't exist")
    }
    # Process the cell barcode whitelist
    if (grepl("\\.gz$", CB_file, ignore.case = TRUE)) {
      con <- gzfile(CB_file, "r")
      sr_barcodes <- readLines(con)
      close(con)
    } else {
      sr_barcodes <- readLines(CB_file)
    }
    # Create cleared barcode file
    cleared_barcodes <- substr(sr_barcodes, 1, 16)
    writeLines(cleared_barcodes, cleared_CB)
    cat("Created cleared barcode file with", length(cleared_barcodes), "barcodes.\n")

  } else{
    # Cleared barcode file already exists
    cat("Cleared barcode file already exists. Using existing file...\n")
  }
  cell_ranger_barcodes <- readLines(cleared_CB)
  # Filter barcodes
  filtered_sg_barcodes <- list()
  for (sgRNA in names(unique_sg_barcodes)) {
    # Extract first 16 bases
    cell_barcodes <- substr(unique_sg_barcodes[[sgRNA]], 1, 16)
    # Keep valid cell barcodes
    valid_indices <- which(cell_barcodes %in% cell_ranger_barcodes)
    valid_barcodes <- cell_barcodes[valid_indices]
    # Store filtered barcodes
    filtered_sg_barcodes[[sgRNA]] <- valid_barcodes
  }

  # Step 4
  cat("Step 4: Extracting unique cell barcodes with high counts...\n")
  for (i in 1:length(gene_list)) {
    gene <- gene_list[i]
    sg <- sg_list[i]
    # count occurences of each barcode
    barcode_counts <- table(filtered_sg_barcodes[[sg]])
    # Filter for barcodes with count > num_counts
    high_count_barcodes <- names(barcode_counts[barcode_counts > num_counts])
    # Write to output file
    output_file <- file.path(outdir, paste0(gene, ".txt"))
    writeLines(high_count_barcodes, output_file)
  }

  cat("Workflow completed successfully!\n")
}
