% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProcessingPipelineFunctions.R
\name{htseq_locusLog2cpm}
\alias{htseq_locusLog2cpm}
\title{Get the log2cpm of a given locus in the htseq output}
\usage{
htseq_locusLog2cpm(
  htseq_filename,
  gene_id,
  markers = c("CNAG_NAT", "CNAG_G418")
)
}
\arguments{
\item{htseq_filename}{filepath to htseq output}

\item{gene_id}{gene id for which to return the log2cpm}

\item{markers}{markers in library. Default CNAG_NAT, CNAG_G418}
}
\value{
log2cpm of a given gene
}
\description{
Get the log2cpm of a given locus in the htseq output
}
