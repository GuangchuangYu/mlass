logisticRegressionResult <- setRefClass("logisticRegrssionResult",
                                  fields=list(
                                  X="matrix",
                                  y="numeric",
                                  theta="matrix"
                                  ),
                                  methods=list(
                                  theta=function(){
                                      .self$theta
                                  },
                                  plot=function(){
                                      xx <- .self$X
                                      yy <- .self$y
                                      tt <- .self$theta

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
                                      return(p)
                                  }
                                  )
                                  )


logisticRegression <- function(X, y, max.iter=20) {
    theta <- rep(0, ncol(X))
    for (i in 1:max.iter) {
        theta <- theta - solve(Hessian(theta, X)) %*% grad(theta,X,y)
    }
    res <- logisticRegressionResult$new(X=X, y=y, theta=theta)
    return(res)
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


