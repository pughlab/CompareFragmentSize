% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format.pvalue.R
\name{format.pvalue}
\alias{format.pvalue}
\title{Format p-value to show on plot}
\usage{
\method{format}{pvalue}(x, digits = 1, lower.value = 1e-08)
}
\arguments{
\item{x}{value to be formatted}

\item{digits}{number of digits to show}

\item{lower.value}{lowest value to show for scientific notation}
}
\value{
expression containing formatted p-value
}
\description{
This function formats a p-value for so that it is easy to read on a plot.
}
\examples{
format.pvalue(0.092);
format.pvalue(0.000006231);
}
