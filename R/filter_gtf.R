#' Filter gtf file with valid transcript_id
#' @param input_file The raw gtf file
#' @param output_file The path to the filtered gtf file
#' @export
filter_gtf <- function(input_file, output_file) {
  fin <- file(input_file, "r")
  out <- file(output_file, "w")
  while (length(line <- readLines(fin, n = 1, warn = FALSE)) > 0) {
    if (grepl("transcript_id", line)) {
      writeLines(line, out)
    }
  }
  close(fin)
  close(out)
  message("Filtered GTF file: ", output_file)
}
