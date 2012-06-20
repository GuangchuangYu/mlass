## logisticRegressionResult <- setRefClass("logisticRegrssionResult",
##                                   fields=list(
##                                   X="matrix",
##                                   y="numeric",
##                                   theta="matrix"
##                                   ),
##                                   methods=list(
##                                   plot=function(){
##                                       xx <- .self$X
##                                       yy <- .self$y
##                                       tt <- .self$theta

##                                       ## compute slope and intercept
##                                       x1 <- c(min(xx[,2]), max(xx[,2]))
##                                       x2 <- -1/tt[3,] * (tt[2,]*x1+tt[1,])
##                                       a <- (x2[2]-x2[1])/(x1[2]-x1[1])
##                                       b <- x2[2]-a*x1[2]


##                                       d <- data.frame(xx, y=factor(yy))
##                                       colnames(d) <- c("V1", "V2", "V3", "label")
##                                       p <- ggplot(d, aes(V2, V3, color=label))+
##                                           geom_point()

##                                       p <- p+geom_abline(slope=a, intercept=b)
##                                       return(p)
##                                   }
##                                   )
##                                   )


## logisticRegression <- function(X, y, max.iter=20) {
##     theta <- rep(0, ncol(X))
##     for (i in 1:max.iter) {
##         theta <- theta - solve(Hessian(theta, X)) %*% grad(theta,X,y)
##     }
##     res <- logisticRegressionResult$new(X=X, y=y, theta=theta)
##     return(res)
## }

## data(ex2data1)
## aa <- logisticRegression(X, y)
## aa$theta
## aa$plot()


##' Class "logisticRegressionResult"
##' This class represents the result of logistic regression.
##'
##'
##' @name logisticRegressionResult-class
##' @aliases logisticRegressionResult-class
##'   "[", logisticRegressionResult-method plot, logisticRegressionResult-method
##'
##' @docType class
##' @slot X x values (a column of 1 was added)
##' @slot y y values
##' @slot theta theta values
##' @seealso \code{\link{logisticRegression}}
##' @keywords classes
##' @author Guangchuang Yu \url{http://ygc.name}
setClass("logisticRegressionResult",
         representation=representation(
         X="matrix",
         y="numeric",
         theta="matrix"
         )
         )
##' Logistic Regression algorithm
##'
##'
##' @title logisticRegression
##' @param X x values
##' @param y labels
##' @param max.iter number of iteration
##' @return A \code{logisticRegressionResult} instance.
##' @importFrom methods new
##' @export
##' @author Yu Guangchuang \url{http://ygc.name}
##' @keywords manip
logisticRegression <- function(X, y, max.iter=15) {
    theta <- rep(0,ncol(X))
    for (i in 1:max.iter) {
        theta <- theta - solve(Hessian(theta, X)) %*% grad(theta,X,y)
    }

    new("logisticRegressionResult",
        X=X,
        y=y,
        theta=theta
        )
}

##' @exportMethod "["
setMethod(
          f="[",
          signature=signature(x="logisticRegressionResult", i="character"),
          function(x, i, j, ...) {
              if (i == "theta"){
                  return(x@theta)
              }
          }
          )


##' plot method for \code{logisticRegressionResult} instance
##'
##' @name plot
##' @docType methods
##' @rdname plot-methods
##' @aliases plot,logisticRegressionResult,ANY-method
##'
##' @title plot method
##' @param object A \code{logisticRegressionResult} instance.
##' @return ggplot object
##' @importFrom stats4 plot
##' @exportMethod plot
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 opts
##' @importFrom ggplot2 geom_abline
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("plot", signature(x="logisticRegressionResult"),
          function(x, title="", xlab="", ylab="") {
              xx <- x@X
              yy <- x@y
              tt <- x@theta

              ## compute slope and intercept
              x1 <- c(min(xx[,2]), max(xx[,2]))
              x2 <- -1/tt[3,] * (tt[2,]*x1+tt[1,])
              a <- (x2[2]-x2[1])/(x1[2]-x1[1])
              b <- x2[2]-a*x1[2]


              d <- data.frame(xx, y=factor(yy))
              colnames(d) <- c("V1", "V2", "V3", "label")
              p <- ggplot(d, aes(V2, V3, color=label))+
                  geom_point()

              p <- p+geom_abline(slope=a, intercept=b)

              p <- p+xlab(xlab)+ylab(ylab)+
                  opts(title=title)
              return(p)
          }
          )


## sigmoid function

g <- function(z) {
    1/(1+exp(-z))
}

## hypothesis function

h <- function(theta, x) {
    g(x %*% theta)
}

## cost function

J <- function(theta, x, y) {
    m <- length(y)
    s <- sapply(1:m, function(i)
                y[i] * log(h(theta, x[i,])) +
                (1-y[i]) * log(1-h(theta, x[i,]))
                )
    j <- -1/m * sum(s)
    return(j)
}


## gradient

grad <- function(theta, x, y) {
    m <- length(y)
    g <- 1/m * t(x) %*% (h(theta,x)-y)
    return(g)
}

## Hessian
Hessian <- function(theta, x) {
    m <- nrow(x)
    H <- 1/m * t(x) %*% x * diag(h(theta, x)) * diag(1-h(theta,x))
    return(H)
}


