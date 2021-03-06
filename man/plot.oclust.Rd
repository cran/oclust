% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rfunctions.R
\name{plot.oclust}
\alias{plot.oclust}
\title{Plots results of the \sQuote{oclust} algorithm.}
\usage{
\method{plot}{oclust}(x, what = c("BIC", "classification", "KL"),
  dimens = NULL, xlab = NULL, ylab = NULL, ylim = NULL,
  addEllipses = TRUE, ...)
}
\arguments{
\item{x}{An \sQuote{oclust} class object obtained by using \code{\link{oclust}}}

\item{what}{A string specifying the type of graph. The options are:
\itemize{
\item \dQuote{BIC}    a plot of BICs for each number of outliers
\item \dQuote{classification}    a plot of the classifications for the optimal solution.
For data with p>2, if more than two \dQuote{dimens} are specified, a pairs plot is produced. If two \dQuote{dimens} are specified,  a coordinate projection plot is produced with the specified \dQuote{dimens}.
Ellipses corresponding to covariances of mixture components are also drawn if \dQuote{addEllipses = TRUE}.
\item \dQuote{KL}    a plot of Kullback-Liebler divergence for each number of outliers}}

\item{dimens}{a vector specifying the dimensions of the coordinate projections}

\item{xlab, ylab}{optional argument specifying axis labels for the classsification plot}

\item{ylim}{optional limits of the y axis of the BIC and KL plots}

\item{addEllipses}{logical indicating whether to include ellipses corresponding
to the covariances of the mixture components}

\item{...}{other graphical parameters}
}
\description{
Plots results of the \sQuote{oclust} algorithm.
}
