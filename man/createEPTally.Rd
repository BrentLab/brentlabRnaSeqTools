% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExperimentSetFunctions.R
\name{createEPTally}
\alias{createEPTally}
\title{Create the EP tally in long form}
\usage{
createEPTally(unfltr_set, qc1_set, iqr_set, tally_conditions)
}
\arguments{
\item{unfltr_set}{the experimental set (a blrs obj), not filtered for any
quality}

\item{qc1_set}{the same experimental set, filtered for qc1 metrics (not iqr)}

\item{iqr_set}{same as above, but with the IQR filter}

\item{tally_conditions}{a base::alist() set of conditions. The members of
the list must be columns in the metadata}
}
\value{
a dataframe in long format detailing the tally of the EP experiment
}
\description{
Create the EP tally in long form
}
