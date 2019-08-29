\name{oclust}
\alias{oclust}
\title{OCLUST Algorithm}
\usage{
oclust(x, o, G, modelNames = NULL, prior = NULL, mc.cores = 1,
  keepAllRes = F, verb = F)
}
\arguments{
\item{x}{A matrix or dataframe with n rows of observations and p columns}

\item{o}{An upperbound for the number of outliers}

\item{G}{The number of clusters}

\item{modelNames}{The model to fit using the Mclust function. Default is NULL (all models).}

\item{prior}{The prior parameter in the Mclust function. Default is NULL.}

\item{mc.cores}{Number of cores to use if running in parallel. Default is 1}

\item{keepAllRes}{A logical value indicating whether to keep the results from all iterations. Default is F.}

\item{verb}{A logical value indicating whether verbose mode is desired, i.e., whether the value of o should be printed as the algorithm proceeds. Default is F.}
}
\value{
oclust returns a list of class oclust with
\item{ data}{The initial data matrix}
\item{ numO }{The predicted number of outliers}
\item{ G }{The number of clusters}
\item{ outs }{The most likely outliers in order of likelihood}
\item{ class}{The classification for the optimal solution}
\item{ pi.g}{The group proportions for the optimal solution}
\item{ mu}{The cluster means for the optimal solution}
\item{ sigma}{The cluster variances for the optimal solution}
\item{KL}{The KL divergence for each iteration, with the first value being for o=0}
\item{BIC}{The BIC for each iteration, with the first value being for o=0}
\item{bic=bic}{The BIC for the optimal solution}
\item{all_results}{(Optional) The parameters for each run if keepAllRes=T. For each, index i+1 corresponds to o=i}
}
\description{
oclust is a trimming method in model-based clustering. It iterates over possible values for the number of outliers and returns the model parameters for the best model as determined by the minimum KL divergence.
}
\references{
Katharine M. Clark and Paul D. McNicholas (2019), Using subset log-likelihoods to trim outliers in Gaussian mixture models. arXiv preprint arXiv:1907.01136.
}
\author{
Katharine M. Clark and Paul D. McNicholas
}
\examples{
\dontshow{
  data(iris)
  iris.o<-oclust(x=iris[1:100,-5],o=3,G=2,modelNames="VVV")
  summary(iris.o)
  plot(iris.o,what="classification")
  plot(iris.o,what="KL")
}
\donttest{
  data(iris)
  iris.o<-oclust(x=iris[,-5],o=10,G=3,modelNames="VVV")
  summary(iris.o)
  plot(iris.o,what="classification")
  plot(iris.o,what="KL")
}
}
\keyword{"multivariate"}
\keyword{"cluster"}
