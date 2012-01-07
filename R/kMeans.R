##' Class "kMeansResult"
##' This class represents the result of kmeans clustering.
##'
##'
##' @name kMeansResult-class
##' @aliases kMeansResult-class
##'   getClusters, kMeansResult-method plot, kMeansResult-method
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


##' Randomly initial centroids for K clusters
##'
##' Randomly selected K features from dataset X
##' @title kMeansInitCentroids
##' @param X data matrix
##' @param K numberic for setting K clusters
##' @return randomly centroids
##' @seealso \code{\link{kMeans}}
##' @author Guangchuang Yu \url{http://ygc.name}
kMeansInitCentroids <- function(X, K) {
    rand.idx <- sample(1:nrow(X), K)
    centroids <- X[rand.idx,]
    return(centroids)
}
##' finding cloest centroids
##'
##'
##' @title findClosestCentroids
##' @param X dataset
##' @param centroids centroids
##' @return cluster index
##' @author Guangchuang Yu \url{http://ygc.name}
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
##' computing centroids
##'
##'
##' @title computeCentroids
##' @param X dataset
##' @param idx cluster index
##' @param K cluster number
##' @return centroids
##' @author Guangchuang Yu \url{http://ygc.name}
computeCentroids <- function(X, idx, K) {
    centroids <- sapply(1:K, function(i)
                        colMeans(X[idx == i,]))
    centroids <- t(centroids)
    return(centroids)
}

##' kMeans algorithm
##'
##' kmeans algorithm
##' @title kMeans
##' @param X dataset
##' @param centers initial centers or K
##' @param max.iter maximum iterations
##' @return A \code{kMeansResult} instance.
##' @importFrom methods new
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
##' @keywords manip
kMeans <- function(X, centers, max.iter = 10){
    if(length(centers) == 1L) {
        K <- centers
        initCentroids <- kMeansInitCentroids(X, K)
    } else {
        initCentroids <- centers
    }

    K <- nrow(initCentroids)
    centroids <- initCentroids
    preCentroids <- centroids

    for (i in 1:max.iter) {
        idx <- findClosestCentroids(X, centroids)
        centroids <- computeCentroids(X, idx, K)
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


##' @exportMethod plot
setGeneric("plot", function(object, ...) standardGeneric("plot"))

##' @exportMethod getClusters
setGeneric("getClusters", function(object) standardGeneric("getClusters"))

##' @exportMethod getCentroids
setGeneric("getCentroids", function(object) standardGeneric("getCentroids"))


##' plot method for \code{kMeansResult} instance
##'
##' @name plot
##' @docType methods
##' @rdname plot-methods
##'
##' @title plot method
##' @param object A \code{kMeans} instance
##' @param trace tracing centroids when algorithm progress
##' @return graph
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("plot", signature(object="kMeansResult"),
          function (object, trace=F, title="", xlab="", ylab="") {
              require(ggplot2)
              X = object@dataset
              colnames(X) <- c("V1", "V2")
              if(ncol(X) != 2) {
                  stop("plot function only support two features in dataset.")
              }
              idx <- object@clusters
              xx <- data.frame(X, cluster=as.factor(idx))
              p <- ggplot(xx, aes(V1, V2))+
                  geom_point(aes(color=cluster))
              if (trace) {
                  preCentroids <- object@traceCentroids
                  preCentroids <- data.frame(preCentroids,
                                             idx=rep(1:3, nrow(preCentroids)/3))
                  p <- p+geom_point(data=preCentroids,
                                    aes(x=V1, y=V2, shape=2)) +
                                        geom_path(data=preCentroids,
                                                  aes(x=V1, y=V2, group=idx)) +
                                                      xlab(xlab) + ylab(ylab) +
                                                          opts(title=title)
              }
              print(p)
          }
          )

##' getClusters method for \code{kMeansResult} instance
##'
##' @name getClusters
##' @docType methods
##' @rdname getClusters-methods
##'
##' @title getCluster method
##' @param object A \code{kMeansResult} instance.
##' @return cluster index
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("getClusters", signature(object="kMeansResult"),
          function(object) {
              clusters <- object@clusters
              return(clusters)
          }
          )

##' getCentroids method for \code{kMeansResult} instance
##'
##' @name getCentroids
##' @docType methods
##' @rdname getCentroids-methods
##'
##' @title getCentroids method
##' @param object A \code{kMeansResult} instance.
##' @return centroids
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("getCentroids", signature(object="kMeansResult"),
          function(object) {
              centroids <- object@centroids
              return(centroids)
          }
          )

