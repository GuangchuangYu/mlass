##' Class "kMeansResult"
##' This class represents the result of kmeans clustering.
##'
##'
##' @name kMeansResult-class
##' @aliases kMeansResult-class
##'   "[",kMeansResult-method
##'   plot,kMeansResult-method
##'
##' @docType class
##' @slot dataset original dataset
##' @slot clusters clustering result
##' @slot centroids final centroids
##' @slot K K clusters
##' @slot traceCentroids tracing the centroids when algorithm progress
##' @seealso \code{\link{kMeans}}
##' @keywords classes
##' @author Guangchuang Yu \url{http://ygc.name}
setClass("kMeansResult",
         representation=representation(
         dataset="matrix",
         clusters="numeric",
         centroids="matrix",
         K="numeric",
         traceCentroids="matrix"
         )
         )


## Randomly initial centroids for K clusters
kMeansInitCentroids <- function(X, K) {
    rand.idx <- sample(1:nrow(X), K)
    centroids <- X[rand.idx,]
    return(centroids)
}


## finding cloest centroids
##' @useDynLib mlass
findClosestCentroids_cpp <- function(X, centroids) {
  idx <- .Call("findClosestCentroids",
              X, centroids,
              package="mlass"
              )
  return(idx)
}


findClosestCentroids <- function(X, centroids) {
    ## finding closest centroids

    # set K
    K <- nrow(centroids)

    idx <- sapply(1:nrow(X), function(i) {
        which.min(
                  sapply(1:K, function(j) {
                      sum(
                          (X[i,]-centroids[j,])^2
                          )
                  })
                  )
    })
    return(idx)
}


## computing centroids
computeCentroids <- function(X, idx, K) {
    centroids <- sapply(1:K, function(i)
                        colMeans(X[idx == i,]))
    centroids <- t(centroids)
    return(centroids)
}

##' @useDynLib mlass
computeCentroids_cpp <- function(X, idx, K) {
  centroids <- .Call("computeCentroids",
                    X, idx, K,
                    package="mlass")
  return(centroids)
}

##' kMeans algorithm
##'
##' kmeans algorithm
##' @title kMeans
##' @param X dataset
##' @param centers initial centers or K
##' @param max.iter maximum iterations
##' @param lang R or CPP
##' @return A \code{kMeansResult} instance.
##' @importFrom methods new
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
##' @keywords manip
kMeans <- function(X, centers, max.iter = 10, lang="CPP"){
    X <- as.matrix(X)
    if(length(centers) == 1L) {
        K <- centers
        initCentroids <- kMeansInitCentroids(X, K)
    } else {
        initCentroids <- as.matrix(centers)
    }

    K <- nrow(initCentroids)
    centroids <- initCentroids
    preCentroids <- centroids

    for (i in 1:max.iter) {

      if (lang == "CPP") {
        idx <- findClosestCentroids_cpp(X, centroids)
        centroids <- computeCentroids_cpp(X, idx, K)
      } else if (lang == "R") {
        idx <- findClosestCentroids(X, centroids)
        centroids <- computeCentroids(X, idx, K)
      }
        preCentroids <- rbind(preCentroids, centroids)
    }

    new("kMeansResult",
        dataset=X,
        clusters=idx,
        centroids=centroids,
        K=K,
        traceCentroids=preCentroids
        )
}


##' plot method for \code{kMeansResult} instance
##'
##' @name plot
##' @docType methods
##' @rdname plot-methods
##' @aliases plot,kMeansResult
##'
##' @title plot method
##' @param x A \code{kMeans} instance
##' @param trace tracing centroids when algorithm progress
##' @return graph
##' @importFrom stats4 plot
##' @exportMethod plot
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_path
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 opts
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("plot", signature(x="kMeansResult"),
          function (x, trace=F, title="", xlab="", ylab="") {
              X = x@dataset

              V1 = colnames(X)[1]
              V2 = colnames(X)[2]
              ## colnames(X) <- c("V1", "V2")

              if(ncol(X) > 2) {
                  warnings("plot function only visualize features in the first two columns.")
              }

              idx <- x@clusters
              xx <- data.frame(X, cluster=as.factor(idx))
              p <- ggplot(xx, aes_string(x=V1, y=V2))+
                  geom_point(aes(color=cluster))
              if (trace) {
                  preCentroids <- x@traceCentroids
                  colnames(preCentroids) <- paste("V", 1:ncol(preCentroids), sep="")
                  K <- x@K
                  preCentroids <- data.frame(preCentroids,
                                             idx=rep(1:K, nrow(preCentroids)/K))
                  p <- p+geom_point(data=preCentroids,
                                    aes(x=V1, y=V2)) +
                                        geom_path(data=preCentroids,
                                                  aes(x=V1, y=V2, group=idx))
                  p <- p + xlab(xlab) + ylab(ylab) +
                      opts(title=title)
              }
              print(p)
          }
          )


##' exportMethod "["
setMethod(
          f="[",
          signature=signature(x="kMeansResult", i="character"),
          definition=function(x,i,j,...) {
              if(i=="clusters")
                  return(x@clusters)
              if(i=="centroids")
                  return(x@centroids)
          }
          )

