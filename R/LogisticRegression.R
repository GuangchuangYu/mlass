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
##' @slot lambda parameter for regularization
##' @seealso \code{\link{logisticRegression}}
##' @keywords classes
##' @author Guangchuang Yu \url{http://ygc.name}
setClass("logisticRegressionResult",
         representation=representation(
         X="matrix",
         y="numeric",
         theta="matrix",
         lambda="numeric"
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
logisticRegression <- function(X, y, max.iter=15, lambda=0) {
    if (lambda == 0) {
        xx <-  X
    } else {
        xx <- apply(X, 1, function(i) mapFeature(i[1], i[2]))
        xx <- t(xx)
    }

    theta <- matrix(rep(0,ncol(xx)), ncol=1)
    for (i in 1:max.iter) {
        theta <- theta - solve(Hessian(theta, xx, lambda)) %*% grad(theta,xx,y, lambda)
    }

    new("logisticRegressionResult",
        X=X,
        y=y,
        theta=theta,
        lambda=lambda
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

              V2 <- V3 <- label <- Var1 <- Var2 <- value <- NULL

              xx <- x@X
              yy <- x@y
              tt <- x@theta
              lambda <- x@lambda

              d <- data.frame(xx, y=factor(yy))
              if (lambda == 0) {
                  colnames(d) <- c("V1", "V2", "V3", "label")
              } else {
                  colnames(d) <- c("V2", "V3", "label")
              }
              p <- ggplot(d, aes(V2, V3, color=label))+
                  geom_point()


              if (lambda == 0) {
                  ## compute slope and intercept
                  x1 <- c(min(xx[,2]), max(xx[,2]))
                  x2 <- -1/tt[3,] * (tt[2,]*x1+tt[1,])
                  a <- (x2[2]-x2[1])/(x1[2]-x1[1])
                  b <- x2[2]-a*x1[2]

                  p <- p+geom_abline(slope=a, intercept=b)
              } else {
                  zz <- logisticPredict(x)
                  p <- p+geom_contour(data=zz, aes(x=Var1,
                                      y=Var2,
                                      z=value),
                                      color="green",
                                      bins=1)
              }
              p <- p+xlab(xlab)+ylab(ylab)+
                  opts(title=title)
              return(p)
          }
          )

##' @importFrom reshape2 melt
logisticPredict <- function(obj){
    ## obj is a "logisticRegressionResult" instance
    x <- obj@X
    theta <- obj@theta
    u <- seq(min(x[,1]), max(x[,1]), len=200)
    v <- seq(min(x[,2]), max(x[,2]), len=200)
    z <- matrix(0, length(u), length(v))
    ## z <- sapply(u, function(i)
    ##             sapply(v, function(j)
    ##                    mapFeature(i,j) %*% theta
    ##                    )
    ##             )
    ## z <- t(z)

    for (i in 1:length(u)) {
        for (j in 1:length(v)) {
            f <- mapFeature(u[i],v[j])
            z[i,j] <- f %*% theta
        }
    }
    rownames(z) <- as.character(u)
    colnames(z) <- as.character(v)
    zz <- melt(z)
    return(zz)
}


## sigmoid function

g <- function(z) {
    1/(1+exp(-z))
}

## hypothesis function

h <- function(theta, x) {
    g(x %*% theta)
}

## cost function

J <- function(theta, x, y, lambda=0) {
    m <- length(y)
    ## s <- sapply(1:m, function(i)
    ##             y[i] * log(h(theta, x[i,])) +
    ##             (1-y[i]) * log(1-h(theta, x[i,]))
    ##             )
    ## j <- -1/m * sum(s)
    j <- -1/m * (
                 y %*% log( h(theta, x) ) +
                 (1-y) %*% log( 1-h(theta, x) )
                 )
    ## regularization
    r <- theta^2
    r[1] <- 0
    j <- j + lambda/(2*m) * sum(r)

    return(j)
}


## gradient

grad <- function(theta, x, y, lambda=0) {
    m <- length(y)
    r <- lambda/m * theta
    r[1] <- 0
    g <- 1/m * t(x) %*% (h(theta,x)-y) + r
    return(g)
}

## Hessian
Hessian <- function(theta, x, lambda=0) {
    m <- nrow(x)
    n <- ncol(x)
    r <- lambda/m * diag(n)
    r[1] <- 0
    H <- 1/m * t(x) %*% x * diag(h(theta, x)) * diag(1-h(theta,x)) + r
    return(H)
}


mapFeature <- function(u, v, degree=6) {
    out <- sapply(0:degree, function(i)
                  sapply(0:i, function(j)
                         u^(i-j) * v^j
                         )
                  )

    out <- unlist(out)
    return(out)
}
