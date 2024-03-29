% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rfunctions.R
\name{findGrossOuts}
\alias{findGrossOuts}
\title{Find Initial Gross Outliers}
\usage{
findGrossOuts(X, minPts = 10, xlim = NULL, elbow = NULL)
}
\arguments{
\item{X}{A data matrix}

\item{minPts}{The minimum number of points in each region of high density. Default is 10}

\item{xlim}{A vector of form c(xmin,xmax) to specify the domain of the plot. Default is NULL, which sets xmax to 10\% of the data size.}

\item{elbow}{An integer specifying the location of the elbow in the plot of Mahalanobis distances. Default is NULL, which returns the plot. If elbow is specified, no plot is produced and the gross outliers are returned.}
}
\value{
findGrossOuts returns a vector with the indices of the gross outliers. One fewer point is returned than the elbow specified.
}
\description{
findGrossOuts uses DBSCAN to find areas of high density. Mahalanobis distance to the closest area of high density is calculated for each point. With no elbow specified, the Mahalonis distances are plotted. If the elbow is specified, the indices of the gross outliers are returned.
}
\details{
The function plots Mahalanobis distance to the closest centre in decreasing order or returns the indices of the gross outliers. The elbow location of the plot provides a good indication as to where the gross outliers end. Running the function first without an elbow specified will plot the Mahalonobis distances. Running it again with the elbow specified will return the outliers. It is recommended to choose the elbow conservatively. If the MDs decrease smoothly, there are no gross outliers. Set elbow=1.
}
