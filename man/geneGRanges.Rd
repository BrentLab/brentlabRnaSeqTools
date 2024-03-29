% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProcessingPipelineFunctions.R
\name{geneGRanges}
\alias{geneGRanges}
\title{Given a GenomicFeatures annotation_db and a gene_id, extract an GRanges
object of the cds}
\usage{
geneGRanges(
  annote_obj_path,
  gene_id,
  id_col = "ID",
  feature_col = "type",
  feature = "cds"
)
}
\arguments{
\item{annote_obj_path}{path to an annotation file parsed by rtracklayer::import}

\item{gene_id}{the ID of a gene in the db. Eg, for cryptococcus CKF44_05222}

\item{id_col}{gene feature column. Default is 'ID'}

\item{feature_col}{feature (col 3) column of annote_obj. Default is 'type'}

\item{feature}{what feature to select. Default is 'cds'}
}
\value{
an IRanges object of a given feature (eg, a gene's cds features)
}
\description{
Given a GenomicFeatures annotation_db and a gene_id, extract an GRanges
object of the cds
}
\references{
rtracklayer::import
}
