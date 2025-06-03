#' @importFrom reticulate import_from_path dict
#' @importFrom basilisk basiliskRun
#' @export
quantify_gene <- function(annotation, outdir, infq, n_process, pipeline = "sc_single_sample", samples = NULL, random_seed = 2024) {
  cat(format(Sys.time(), "%X %a %b %d %Y"), "quantify genes \n")

  if (grepl("\\.gff3?(\\.gz)?$", annotation)) {
    warning("Annotation in GFF format may cause errors. Please consider using GTF formats.\n")
  }

  genome_bam <- list.files(outdir)[grepl("_?align2genome\\.bam$", list.files(outdir))]
  cat("Found genome alignment file(s): ")
  cat(paste0("\t", paste(genome_bam, collapse = "\n\t"), "\n"))

  if (length(genome_bam) != 1 && grepl("single_sample", pipeline)) {
    stop("Incorrect number of genome alignment files found.\n")
  }

  tryCatch(
    {
      basiliskRun(
        env = FLAMES:::flames_env, fun = function(annotation, outdir, pipeline, n_process, infq, samples, random_seed) {
          python_path <- system.file("python", package = "FLAMES")
          count <- reticulate::import_from_path("count_gene", python_path)
          count$quantification(annotation, outdir, pipeline, n_process, infq = infq, sample_names = samples, random_seed = random_seed)
        },
        annotation = annotation,
        outdir = outdir,
        pipeline = pipeline,
        n_process = n_process,
        infq = infq,
        samples = samples,
        random_seed = random_seed
      )
    },
    error = function(e) {
      # Capture the Python error using py_last_error()
      py_error <- reticulate::py_last_error()
      if (!is.null(py_error)) {
        py_error_message <- py_error$message
        # Print the actual function call
        cat(annotation, outdir, pipeline, n_process, infq, samples, random_seed)
        stop("Error when quantifying genes:\n", py_error_message)
      } else {
        stop("Error when quantifying genes:\n", e$message)
      }
    }
  )
}
