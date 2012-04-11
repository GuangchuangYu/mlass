\name{svmTrain}
\alias{svmTrain}
\title{svmTrain}
\usage{
  svmTrain(X, Y, C, kernelFunction = "linearKernel",
    tol = 1e-04, max.iter = 20, verbose = FALSE)
}
\arguments{
  \item{X}{dataset}

  \item{Y}{label}

  \item{C}{control penalty}

  \item{kernelFunction}{kernel function}

  \item{tol}{tolerance}

  \item{max.iter}{max iteration}

  \item{verbose}{logical parameter}
}
\value{
  A \code{svmResult} instance
}
\description{
  SVM training
}
\details{
  SMO algorithm
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}
\keyword{manip}

