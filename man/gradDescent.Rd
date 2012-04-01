\name{gradDescent}
\alias{gradDescent}
\title{gradDescent}
\usage{
  gradDescent(X, y, theta, alpha = 0.01, max.iter = 1500)
}
\arguments{
  \item{X}{x values (a column of 1 was added)}

  \item{y}{y values}

  \item{theta}{initial theta values}

  \item{alpha}{learning rate}

  \item{max.iter}{number of iteration}
}
\value{
  A \code{linearRegressionResult} instance.
}
\description{
  Gradient descent algorithm for linear regression
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}
\keyword{manip}

