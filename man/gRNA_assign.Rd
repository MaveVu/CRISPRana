% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gRNA_assign.R
\name{gRNA_assign}
\alias{gRNA_assign}
\title{Generaate gRNA files containing gRNA-cell barcode mapping information}
\usage{
gRNA_assign(
  ref_sgRNA,
  gRNA_R1,
  gRNA_R2,
  CB_file,
  cleared_CB,
  num_counts = 10,
  outdir
)
}
\arguments{
\item{ref_sgRNA}{The reference file containing gRNA names and sequences}

\item{gRNA_R1}{The Read 1 data containing cell barcodes and UMIs}

\item{gRNA_R2}{The Read 2 data containing gRNA sequences}

\item{CB_file}{The cell barcode whitelist file}

\item{cleared_CB}{The file containing cleared cell barcodes (no "-1")}

\item{num_counts}{The threshold for filtering low-counts cell barcodes}

\item{outdir}{The directory containing gRNA files}
}
\value{
A directory containing gRNA files
}
\description{
Generaate gRNA files containing gRNA-cell barcode mapping information
}
