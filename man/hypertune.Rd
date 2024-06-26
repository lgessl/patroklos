% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tune.R
\name{hypertune}
\alias{hypertune}
\title{Build a fitter with integrated cross-validation from a fitter}
\usage{
hypertune(fitter, select = FALSE)
}
\arguments{
\item{fitter}{A patroklos-compliant fitter.}

\item{select}{logical. If \code{TRUE}, fitter will only return the model that
minimizes the error. If \code{FALSE}, fitter will return a \code{ptk_hypertune} object
with all models and their errors.}
}
\value{
A patroklos-compliant fitter with integrated cross-validation tuning.
}
\description{
This function operator takes a patroklos-compliant fitter and
builds a patroklos-compliant fitter with integrated cross validation from
it.
}
