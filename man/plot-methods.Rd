\docType{methods}
\name{plot}
\alias{plot}
\alias{plot,kMeansResult}
\alias{plot,linearRegressionResult,ANY-method}
\title{plot method}
\arguments{
  \item{x}{A \code{kMeans} instance}

  \item{trace}{tracing centroids when algorithm progress}

  \item{object}{A \code{linearRegressionResult} instance.}

  \item{x}{A \code{svmResult} instance}

  \item{X}{data matrix}

  \item{y}{class label, 0 or 1}

  \item{type}{one of linear or nonlinear}

  \item{title}{title}

  \item{xlab}{xlab}

  \item{ylab}{ylab}
}
\value{
  graph

  ggplot object

  ggplot2 graph object
}
\description{
  plot method for \code{kMeansResult} instance

  plot method for \code{linearRegressionResult} instance

  plot method for \code{svmResult} instance
}
\author{
  Guangchuang Yu \url{http://ygc.name}

  Guangchuang Yu \url{http://ygc.name}

  Guangchuang Yu \url{http://ygc.name}
}

