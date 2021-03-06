% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bundle_theta_y_sample.R
\name{bundle_theta_y_sample}
\alias{bundle_theta_y_sample}
\title{Bundle result of Metropolis sampling into a single convenient list}
\usage{
bundle_theta_y_sample(
  theta_yList,
  logLikVect,
  x,
  Y,
  hp,
  varNames = NA,
  known = NA
)
}
\arguments{
\item{theta_yList}{A list of samples of theta_y (samples are a list, not a vector)}

\item{x}{Vector of indepenent variable observations}

\item{Y}{Matrix of dependent variable observations}

\item{hp}{Hyperparameters}

\item{varNames}{Vector of variable names [optional]}

\item{known}{List object of known simulation values [optional]}

\item{logLikeVect}{A vector of likelihoods calculated for each sample}
}
\value{
The bundled list
}
\description{
\code{bundles_theta_y_sample} bundles theta_yList, x, Y, hp, varNames, and known into a single convenient list. This bundled result is used for displaying and saving results.
}
\seealso{
simConstHaz, fitConstHaz, fitLinHaz
}
\author{
Michael Holton Price <MichaelHoltonPrice@gmail.com>
}
