% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_initialization_h.R
\name{get_initialization}
\alias{get_initialization}
\title{get_initialization}
\usage{
get_initialization(input_data, purity, phi = c(), kappa = 5, alpha = 0.05)
}
\arguments{
\item{input_data}{list: List of 7: $S: int, $N: int, $karyotype: num (0 or 1), $seg_assignment: num, $peaks:List of N of num (1:2), $NV: num, $DP: num}

\item{purity}{num: (type = double) sample purity between 0 and 1}

\item{phi}{parameters of dirichlet reparametrization uniformly initialized}

\item{kappa}{parameters of dirichlet reparametrization}

\item{alpha}{num: (type double) confidence interval level to choose the data that fall in the expected binomial intervals}
}
\value{
inits list: List of 4: $w: num (1:S, 1:3), $tau: num (1:K), $phi: num (1:K), $kappa: num
}
\description{
Perform c-means on the proportions of mutations that fall in to the first or second binomial interval (in accepted mutations) and retrieve the centroids and u matrix to be used as initialization parameters.
}
\keyword{cluster}
