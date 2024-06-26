% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_timing.R
\name{plot_timing}
\alias{plot_timing}
\title{Plot timing of clonal peaks in cancer genome sequencing data}
\usage{
plot_timing(fit_results, segments, colour_by = "karyotype", ref = "GRCh38")
}
\arguments{
\item{fit_results}{A list containing the results of the \code{fit_timing} function, specifically \code{summarized_results}.}

\item{segments}{A data frame containing segment information with columns \code{chr}, \code{from}, \code{to}, \code{Major}, and \code{minor}.}

\item{colour_by}{A character string specifying the variable to color the plot by. Default is "karyotype".}

\item{ref}{Reference genome desired. Either 'GRCh38' or 'hg19'}
}
\value{
A ggplot object showing the timing of clonal peaks.
}
\description{
This function generates a plot showing the timing of clonal peaks in cancer genome sequencing data.
}
