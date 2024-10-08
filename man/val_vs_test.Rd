% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/val_vs_test.R
\name{val_vs_test}
\alias{val_vs_test}
\title{Compare validation and test error}
\usage{
val_vs_test(
  model_list,
  data,
  error_fun,
  regex1 = NULL,
  regex2 = NULL,
  name1 = NULL,
  name2 = NULL,
  legendtitle1 = "spot 1",
  legendtitle2 = "spot 2",
  correlation_label = TRUE,
  file = NULL,
  return_type = c("ggplot", "tibble"),
  plot_theme = ggplot2::theme_minimal(),
  colors = NULL,
  width = 7,
  height = 4,
  quiet = FALSE
)
}
\arguments{
\item{model_list}{list of \code{Model} objects.}

\item{data}{\code{Data} object. Its \code{cohort} attribute must be set. Calculate the test error on the
\code{cohort} of \code{data}.}

\item{error_fun}{function. Error function to calculate the error between the true and predicted
response. For its interface, see \code{val_error_fun_prototype()}.}

\item{regex1}{character vector. Regular-expression patterns. For every model, we assign it the
group in \code{name1} corresponding to the first pattern its name matches and color it accordingly.}

\item{regex2}{character vector. Regular-expression patterns. For every model, we assign it the
group in \code{name2} corresponding to the first pattern its name matches and shape it accordingly.}

\item{name1, name2}{character vectors. Names of the groups as defined via \code{regex1} and \code{regex2}.
Therefore, \code{length(name1) == length(regex1)} and \code{length(name2) == length(regex2)}. If \code{name1}
or \code{name2} is \code{NULL}, we use the corresponding \code{regex1} or \code{regex2} as the name.}

\item{legendtitle1, legendtitle2}{string. Legend titles for the grouping according to \code{regex1}
and \code{regex2}.}

\item{correlation_label}{logical. Whether to show the correlation coefficient between the
validation and test error as a label in the plot.}

\item{file}{string. File name to save the plot to.}

\item{return_type}{string. Either "ggplot" or "tibble". See return section for details.}

\item{plot_theme}{ggplot2 theme. Theme to apply to the plot. If it's \code{NULL}, use the ggplot2
default theme.}

\item{colors}{character vector. Colors used for points.}

\item{width, height}{numeric. Width and height of the stored plot in inches.}

\item{quiet}{logical. Whether to suppress messages.}
}
\value{
ggplot object if \code{return_type == "ggplot"} or the tibble underlying the plot if
\code{return_type == "tibble"}.
}
\description{
Plot validation error against test error for a list of models.
}
