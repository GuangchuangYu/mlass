##' Class "linearRegressionResult"
##' This class represents the result of linear regression.
##'
##'
##' @name linearRegressionResult-class
##' @aliases linearRegressionResult-class
##'  "[", lineareRegressionResult-method plot, linearRegressionResult-method
##'
##' @docType class
##' @slot X x values (a column of 1 was added)
##' @slot y y values
##' @slot theta theta values
##' @seealso \code{\link{gradDescent}}
##' @keywords classes
##' @author Guangchuang Yu \url{http://ygc.name}
setClass("linearRegressionResult",
         representation=representation(
         X="matrix",
         y="numeric",
         theta="matrix"
         )
         )




##' computing the cost J
##'
##' The gradient descent algorithm was to minimize the cost function.
##' @title computeCost
##' @param X x values (a column of 1 was added)
##' @param y y values
##' @param theta theta values
##' @return cost J
##' @author Guangchuang Yu \url{http://ygc.name}
computeCost <- function(X, y, theta, lambda=0) {
    ## number of training example.
    m <- length(y)

    r <- theta^2
    r[1] <- 0
    h <- theta %*% t(X)
    hh <- t(h) -y
    J <- 1/(2*m) * sum(hh^2) + lambda*sum(r)

    return(J)
}

##' linear regression
##'
##'
##' @title linearRegression
##' @param X x values (a column of 1 was added)
##' @param y y values
##' @param theta initial theta values
##' @param alpha learning rate
##' @param max.iter number of iteration
##' @return A \code{linearRegressionResult} instance.
##' @importFrom methods new
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
##' @keywords manip
linearRegression <- function(X,y, theta, alpha=0.01, lambda=0, max.iter=1000, normalEqn=FALSE, featureNormalize=FALSE) {
    if (normalEqn) {
        theta <- normalEqn(X, y)
    } else {
        if (featureNormalize) {
            xx <- featureNormalize(X)
        } else {
            xx <- X
        }
        theta <- gradDescent(xx, y, theta, alpha, max.iter, lambda)
    }
    new("linearRegressionResult",
        X=X,
        y=y,
        theta=theta
        )
}

## data(ex1data3)
x <- mapFeature(X)
theta <- rep(0, ncol(x))
aa <- linearRegression(x, y, theta, lambda=0, alpha=0.01)
x.test <- seq(-1,1, 0.001)
y.test <- mapFeature(x.test) %*% t(aa["theta"])
d <- data.frame(x=X, y=y)
p <- ggplot(data=d, aes(x=x,y=y))+geom_point()
dd <- data.frame(x.test=x.test, y.test=y.test)
p+geom_line(data=dd, aes(x=x.test, y=y.test))


## Normal equation
normalEqn <- function(X, y, lambda=0) {
    n <- ncol(X)
    ## extra regularization terms
    r <- lambda * diag(n)
    r[1,1] <- 0
    theta <- solve(t(X) %*% X + r) %*% t(x) %*% y
    return(theta)
}

## Gradient descent algorithm
gradDescent <- function(X, y, theta, alpha, max.iter, lambda=0) {
    ## number of training examples
    m <- length(y)

    for (i in 1:max.iter) {
        tt <- theta
        tt[1] <- 0
        h <- theta %*% t(X)
        ## h <- X %*% t(theta)
        ## dj <- t(h-y) %*% X + lambda * tt
        dj <- t(t(h)-y) %*% X + lambda * tt ## derivative of cost J.
        theta <- theta - alpha * dj/m
        ## print(theta)
    }
    return(theta)
}

gradDescent2 <- function( x, y,theta, alpha=0.1, niter=1000, lambda=1) {
    m <- length(y)
    for (i in 1:niter) {
        tt <- theta
        tt[1] <- 0
        dj <- 1/m * (t(h(theta,x)-y) %*% x + lambda * tt)
        theta <- theta - alpha * dj

    }

    return(theta)
}


## @exportMethod getTheta
## setGeneric("getTheta", function(object) standardGeneric("getTheta"))

## getTheta method for \code{linearRegressionResult} instance
##
## @name getTheta
## @docType methods
## @rdname getTheta-methods
##
## @title getTheta method
## @param object A \code{linearRegressionResult} instance.
## @return theta
## @author Guangchuang Yu \url{http://ygc.name}
## setMethod(getTheta, signature(object="linearRegressionResult"),
##          function (object) {
##              theta <- object@theta
##              return(theta)
##          }
##          )


##' @exportMethod "["
setMethod(
          f="[",
          signature=signature(x="linearRegressionResult", i="character"),
          function(x, i,j,...) {
              if (i == "theta") {
                  return(x@theta)
              }

          }
          )


##' plot method for \code{linearRegressionResult} instance
##'
##' @name plot
##' @docType methods
##' @rdname plot-methods
##' @aliases plot,linearRegressionResult,ANY-method
##'
##' @title plot method
##' @param object A \code{linearRegressionResult} instance.
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
setMethod("plot",signature(x="linearRegressionResult"),
          function(x, title="", xlab="", ylab="") {
              X <- x@X
              y <- x@y
              theta <- x@theta
              pg <- ggplot()+ aes(X[,2],y) +
                  geom_point()

              ## predicted <- as.vector(theta %*% t(X))
              ## predicted.df <- data.frame(x=X[,2], y=predicted)
              ## pg <- pg+geom_line(data=predicted.df,
              ##                   aes(x=x,y=y, color="red")) +
              ##                       opts(legend.position="none")
              pg <- pg +
                  geom_abline(intercept=theta[1],
                              slope=theta[2],
                              colour="red")

              pg <- pg + xlab(xlab) + ylab(ylab) +
                  opts(title=title)
              return(pg)
          }
          )

## normailze features using Z-score
featureNormalize <- function(x) {
    for (i in 2:ncol(x)) {
        x[,i] <- (x[,i] - mean(x[,i]))/sd(x[,i])
    }
    return(x)
}


mapFeature <- function(x, degree=5){
    return(sapply(0:degree, function(i) x^i))
}
