% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tune.R
\name{multitune}
\alias{multitune}
\title{Tune multiple hyperparameters with a single call to a fitter}
\usage{
multitune(fitter, select = FALSE)
}
\arguments{
\item{fitter}{A fitting function that fulfills the \code{\link[=fitter_prototype]{fitter_prototype()}} interface.}

\item{select}{logical. If \code{TRUE}, the returned \code{multitune_obj} S3 object will only contain the
validated predictions and the parameters of the picked model among the fitted models.}
}
\value{
A fitter following the interface of \code{\link[=fitter_prototype]{fitter_prototype()}} documented in
\code{\link[=multitune_output_prototype]{multitune_output_prototype()}}.
}
\description{
Decorate a fitter that can only tune one combination of hyperparameters at a time
into a fitter that can tune multiple combinations at once.
}
