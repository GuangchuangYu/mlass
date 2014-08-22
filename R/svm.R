## support vector machine

## SVM classifier using a simplified version of sequential minimal optimization (SMO) algorithm.

## This is a simplified version of the SMO algorithm for training SVMs.


## We’d like to find the point closest to the separating hyperplane
## and make sure this is as far away from the separating line as possible.
## This is known as margin. We want to have the greatest possible margin,
## because if we made a mistake or trained our classifier on limited data,
## we’d want it to be as robust as possible.



##' Class "svmResult"
##' This class represents the result of svm training model.
##'
##'
##' @name svmResult-class
##' @aliases svmResult-class
##'     "[",svmResult-method
##'     plot,svmResult-method
##'
##' @docType class
##' @slot X dataset
##' @slot y label
##' @slot kernelFunction kernnelFunction
##' @slot w w
##' @slot b b
##' @slot alphas alphas
##' @seealso \code{\link{svmTrain}}
##' @keywords classes
##' @author Guangchuang Yu \url{http://ygc.name}
setClass("svmResult",
         representation=representation(
         X="matrix",
         y="numeric",
         kernelFunction="character",
         w="matrix",
         b="numeric",
         alphas="numeric"
         )
         )




##' SVM training
##'
##' SMO algorithm
##' @title svmTrain
##' @param X dataset
##' @param Y label
##' @param C control penalty
##' @param kernelFunction kernel function
##' @param tol tolerance
##' @param max.iter max iteration
##' @param verbose logical parameter
##' @return A \code{svmResult} instance
##' @importFrom methods new
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
##' @keywords manip
svmTrain <- function(X, Y, C, kernelFunction="linearKernel",
                     tol=1e-4, max.iter=20, verbose=FALSE) {
    #####################################################################
    ##                                                                 ##
    ##   X is the matrix of training examples,                         ##
    ##     each row is a training example,                             ##
    ##     and jth column holds the jth feature.                       ##
    ##                                                                 ##
    ##   Y is a vector containing 1 for positive examples              ##
    ##     and 0 for negative examples.                                ##
    ##                                                                 ##
    ##   C is the standard SVM regularization parameter.               ##
    ##                                                                 ##
    ##   tol is the tolerance value.                                   ##
    ##                                                                 ##
    ##   max.iter controls the number of iterations over the dataset   ##
    ##     (without changes to alpha) before the algorithm quits.      ##
    ##                                                                 ##
    #####################################################################


    if (verbose) {
        print("Training...")
        pb <- txtProgressBar(min = 0, max = max.iter, style = 2)
    }

    ## data parameter
    m <- nrow(X)
    #n <- ncol(X)

    ## labels
    ## Map 0 to -1
    Y[Y==0] <- -1

    # variables
    alphas <- rep(0, m)
    b <- 0
    E <- rep(0, m)
    iter <- eta <- L <- H <- 0

    ## transformations of original data to map into new space
    if(kernelFunction == "linearKernel") {
        K <- X %*% t(X)
    } else if (kernelFunction == "gaussianKernel") {
        K <- matrix(0, ncol=m, nrow=m)
        for (i in 1:m) {
            for (j in i:m) {
                K[i,j] <- gaussianKernel(X[i,], X[j,])
                K[j,i] <- K[i,j]
            }
        }
    }

    while (iter < max.iter) {
        num_changed_alphas <- 0
        for (i in 1:m){
            E[i] <- b + sum(alphas * Y * K[,i]) - Y[i]
            if( (Y[i]*E[i] < -tol & alphas[i] < C) || (Y[i] >  tol & alphas[i] > 0) ) {
                ## if the error E[i] is large
                ## the alpha corresponding to this data instance can be optimized

                j <- ceiling(m * runif(1))
                while (j == i) { # make sure i != j
                    j <- ceiling(m * runif(1))
                }

                E[j] <- b + sum(alphas * Y * K[,j]) - Y[j]

                ## save old alphas
                alpha.i.old <- alphas[i]
                alpha.j.old <- alphas[j]

                if (Y[i] == Y[j]) {
                    L <- max(0, alphas[j] + alphas[i] -C)
                    H <- min(C, alphas[j] + alphas[i])
                } else {
                    L <- max(0, alphas[j] - alphas[i])
                    H <- min(C, C + alphas[j] - alphas[i])
                }

                if (L == H) {
                    ## continue to next i
                    next
                }

                ## compute eta
                eta <- 2 * K[i,j] - K[i,i] - K[j,j]
                if (eta >= 0) {
                    ## continue to next i
                    next
                }

                ## compute and clip new value for alpha j
                alphas[j] <- alphas[j] - Y[j] * (E[i] - E[j])/eta
                ## clip
                alphas[j] = min(H, alphas[j])
                alphas[j] = max(L, alphas[j])

                ## check if change in alpha is significant
                if(abs(alphas[j] - alpha.j.old) < tol) {
                    alphas[j] <- alpha.j.old
                    next
                }

                ## determine value for alpha i
                alphas[i] <- alphas[i] + Y[i] * Y[j] * (alpha.j.old - alphas[j])

                ## compute b1 and b2
                b1 <- b - E[i] +
                    - Y[i] * (alphas[i] - alpha.i.old) * K[i,j] +
                        - Y[j] * (alphas[j] - alpha.j.old) * K[i,j]
                b2 <- b - E[j] +
                    - Y[i] * (alphas[i] - alpha.i.old) * K[i,j] +
                        - Y[j] * (alphas[j] - alpha.j.old) * K[i,j]

                ## compute b
                if ( alphas[i] > 0 & alphas[i] < C) {
                    b <- b1
                } else if (alphas[j] > 0 & alphas[j] < C) {
                    b <- b2
                } else {
                    b <- (b1+b2)/2
                }

                num_changed_alphas <- num_changed_alphas + 1
            }
        }

        if (num_changed_alphas == 0) {
            iter <- iter + 1
        } else {
            iter <-  0
        }
        if (verbose) {
            setTxtProgressBar(pb, iter)
        }
    }

    if(verbose) {
        close(pb)
        print("done")
    }

    idx <- alphas > 0
    new("svmResult",
        X=X[idx,],
        y=Y[idx],
        kernelFunction=kernelFunction,
        b=b,
        alphas=alphas[idx],
        w=t(alphas * Y) %*% X
        )
}

##' plot method for \code{svmResult} instance
##'
##' @name plot
##' @docType methods
##' @rdname plot-methods
##' @aliases plot,svmResult
##'
##' @title plot method
##' @param x A \code{svmResult} instance
##' @param X data matrix
##' @param y class label, 0 or 1
##' @param type one of linear or nonlinear
##' @param title title
##' @param xlab xlab
##' @param ylab ylab
##' @return ggplot2 graph object
##' @importFrom stats4 plot
##' @exportMethod plot
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 geom_contour
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggtitle
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("plot", signature(x="svmResult"),
          function(x, X, y, type, title="", xlab="", ylab=""){
              model <- x
              match.arg(type, c("linear", "nonlinear"))
              d <- data.frame(X, y=as.factor(y))
              V <- colnames(X)
              p <- ggplot(data=d, aes_string(x=V[1], y=V[2])) +
                  geom_point(aes(colour=y))
              w <- model["w"]
              b <- model["b"]
              if(type == "linear") {
                  p <- p+geom_abline(intercept = -b/w[2],
                                     slope = -w[1]/w[2],
                                     colour = "red")
              } else {
                  xr <- seq(min(X[,1]), max(X[,1]), length.out=100)
                  yr <- seq(min(X[,2]), max(X[,2]), length.out=100)
                  mg <- meshgrid(xr, yr)
                  X1 <- mg[[1]]
                  X2 <- mg[[2]]

                  vals <- matrix(0, ncol=ncol(X1), nrow=nrow(X1))

                  for (i in 1:ncol(X1)){
                      thisX <- cbind(X1[,i], X2[,i])
                      vals[,i] <- svmPredict(model, thisX)
                  }

                  vm <- data.frame(x=as.vector(X1),
                                   y= as.vector(X2),
                                   z=as.vector(vals))

                  p <- p+geom_contour(data=vm, aes(x=x,y=y,z=z))
              }

              p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)

              print(p)
          }
          )

##' exportMethod "["
setMethod(
          f="[",
          signature=signature(x="svmResult", i="character"),
          definition=function(x,i,j...) {
              if(i=="w")
                  return(x@w)
              if(i=="b")
                  return(x@b)
              if(i=="alphas")
                  return(x@alphas)
              if(i=="kernelFunction")
                  return(x@kernelFunction)
              if(i=="X")
                  return(x@X)
              if(i=="y")
                  return(x@y)
          }
          )

gaussianKernel <- function(x1,x2, sigma=0.1) {
    ## Gaussian kernel is a similarity function that
    ## measures the "distance" between a pair of examples.
    sim <- exp(-sum((x1-x2)^2)/(2*sigma^2))
    return(sim)
}

svmPredict <- function(model, X) {
    ## returns a vector of predictions using a trained SVM model (svmTrain).
    ## X is a m x n matrix where there each example is a row.
    ## model is a svm model returned from svmTrain.

    ## output pred is a vector of length m of predictions of {0, 1} values.

    if (ncol(X) == 1) {
        ## examples should be in rows
        X <- t(X)
    }

    m <- nrow(X)
    p <- rep(0, m)
    pred <- rep(0, m)
    kernelFunction <- model["kernelFunction"]
    if (kernelFunction == "linearKernel") {
        p <- X %*% model["w"] + model["b"]
    } else if (kernelFunction == "gaussianKernel") {

        ##    X1 <- rowSums(X^2)
        ##    X2 <- rowSums(model$X^2)
        ##    K <- X1 + X2 - 2* X %*% t(model$X)
        ##    K <- gaussianKernel(1,0) ^K
        ##    K <- model$y * K
        ##    K <- model$alphas * K
        ##    p <- rowSums(K)

        alphas <- model["alphas"]
        for (i in 1:m) {
            prediction <- 0
            for (j in 1:nrow(model["X"])) {
                prediction <- prediction +
                    alphas[j] * model["y"][j] *
                        gaussianKernel(X[i,], model["X"][j,])
            }
            p[i] <- prediction + model["b"]
        }

    }
    pred[p >= 0] <- 1
    pred[p < 0] <- 0
    return(pred)
}


## equivalent with meshgrid in octave
meshgrid <- function(a, b) {
    list(
         x <- outer(b*0, a, FUN="+"),
         y <- outer(b, a*0, FUN="+")
         )
}

