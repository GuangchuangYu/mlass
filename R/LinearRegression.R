##' Class "linearRegressionResult"
##' This class represents the result of linear regression.
##'
##'
##' @name linearRegressionResult-class
##' @aliases linearRegressionResult-class
##'   getTheta, lineareRegressionResult-method plot, linearRegressionResult-method
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
computeCost <- function(X, y, theta) {
    ## number of training example.
    m <- length(y)

    h <- theta %*% t(X)
    hh <- t(h) -y
    J <- 1/(2*m) * sum(hh^2)

    return(J)
}

##' Gradient descent algorithm for linear regression
##'
##'
##' @title gradDescent
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
gradDescent <- function(X,y, theta, alpha=0.01, max.iter=1500) {
    ## number of training examples
    m <- length(y)

    for (i in 1:max.iter) {
        h <- theta %*% t(X)
        dj <- t(t(h)-y) %*% X ## derivative of cost J.
        theta <- theta - alpha * dj/m
    }

    new("linearRegressionResult",
        X=X,
        y=y,
        theta=theta
        )
}

##' @exportMethod getTheta
setGeneric("getTheta", function(object) standardGeneric("getTheta"))

##' getTheta method for \code{linearRegressionResult} instance
##'
##' @name getTheta
##' @docType methods
##' @rdname getTheta-methods
##'
##' @title getTheta method
##' @param object A \code{linearRegressionResult} instance.
##' @return theta
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod(getTheta, signature(object="linearRegressionResult"),
          function (object) {
              theta <- object@theta
              return(theta)
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
              require(ggplot2)
              X <- x@X
              y <- x@y
              theta <- x@theta
              pg <- ggplot()+ aes(X[,2],y) +
                  geom_point() +
                      xlab(xlab)+ylab(ylab) +
                          opts(title=title)
              ##predicted <- as.vector(theta %*% t(X))
              ##predicted.df <- data.frame(x=X[,2], y=predicted)
              ##pg <- pg+geom_line(data=predicted.df,
              ##                   aes(x=x,y=y, color="red")) +
              ##                       opts(legend.position="none")
              pg <- pg +
                  geom_abline(intercept=theta[1],
                              slope=theta[2],
                              colour="red")
              return(pg)
          }
          )
