% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/training_camp.R
\name{training_camp}
\alias{training_camp}
\title{A wrapper around \code{Model$fit()}}
\usage{
training_camp(model_list, data, skip_on_error = TRUE, quiet = FALSE)
}
\arguments{
\item{model_list}{list of Model objects. The models to fit.}

\item{data}{Data object. The data set to fit on.}

\item{skip_on_error}{logical. Whether to skip to the next model if an error
occurs while fitting a model.}

\item{quiet}{logical. Whether to suppress messages. Default is \code{FALSE}.}
}
\description{
Given one data set and a list of \code{Model}s, call \code{Model$fit()} for
all of them. If an error occurs while fitting a model, you can skip to the next
model.
}
