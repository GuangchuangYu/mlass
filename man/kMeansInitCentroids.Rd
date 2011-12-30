\name{kMeansInitCentroids}
\alias{kMeansInitCentroids}
\title{kMeansInitCentroids}
\usage{
kMeansInitCentroids(X, K)
}
\arguments{
  \item{X}{data matrix}

  \item{K}{numberic for setting K clusters}
}
\value{
  randomly centroids
}
\description{
  Randomly initial centroids for K clusters
}
\details{
  Randomly selected K features from dataset X
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}
\seealso{
  \code{\link{kMeans}}
}
