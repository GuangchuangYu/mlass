distance <- function(x,y, method="euclidean") {
    switch(method,
           "euclidean" = sqrt(sum((x-y)^2))
           ## implement other methods here.
           )
}

getDistance <- function(mat, method="euclidean") {
    nr <- nrow(mat)
    d <- matrix(NA, ncol=nr, nrow=nr)
    colnames(d) <- rownames(mat)
    rownames(d) <- rownames(mat)

    for (i in 1:nr) {
        for (j in i:nr) {
            d[i,j] <- distance(mat[i,], mat[j,], method)
            if (i != j)
                d[j,i] <- d[i,j]
        }
    }
    return(d)
}

getClosest <- function(d) {
    nr <- nrow(d)
    h <- min(d[lower.tri(d)])
    node <- names(unlist(sapply(1:nr, function(i)
                                which(d[i,] == h))))
    hc <- list(node=node, h=h)
    return(hc)
}

linkage <- function(d, method, cn) {
    sc <- getClosest(d)
    fun <- switch(method,
                  "complete" = max,
                  "single"   = min
                  )
    node <- sc$node
    h <- apply(d[node,], 2, fun)
    h <- h[! names(h) %in% node]
    d <- d[!rownames(d) %in% node,]
    if (is.null(dim(d))) {
        nm <- names(d)
        d <- matrix(d, nrow=1)
        colnames(d) <- nm
        rownames(d) <- nm[!nm %in% node]
    }
    rn <- rownames(d)
    d <- d[,!colnames(d) %in% node]
    if(is.null(dim(d))) {
        d <- matrix(d)
        rownames(d) <- rn
        colnames(d) <- rn
    }
    d <- rbind(d, h)
    d <- cbind(d, h=c(h,0))
    rownames(d)[rownames(d) == "h"] <- cn
    colnames(d)[colnames(d) == "h"] <- cn
    result <- list(subCluster=sc, d=d)
    return(result)
}

hcluster <- function(mat,
                     method="complete",
                     dist.method="euclidean") {
    d <- getDistance(mat, dist.method)
    nr <- nrow(d)
    cls <- list()
    for (i in 1:(nr-1)) {
        cn <- paste("c", i, sep="_")
        sc <- linkage(d, method, cn)
        cls[[i]] <- sc$subCluster
        names(cls)[i] <- cn
        d <- sc$d
    }
    result <- list(clusters=cls,
                   data=mat,
                   method=method,
                   dist.method=dist.method)
    return(result)
}

plotting_hcluster <- function(hclusterResult, main="Cluster Dendrogram", xlab="", ylab="Height") {
    cls <- hclusterResult$clusters
    labels <- rownames(hclusterResult$data)
    tn <- sapply(cls, function(i) rev(i$node))
    h <- sapply(cls, function(i) i$h)
    nl <- tn[,ncol(tn)]
    idx <- which(! nl %in% labels)
    while(length(idx)) {
        i <- idx[1]
        if (i > 1) {
            start <- nl[1:(i-1)]
        } else {
            start <- c()
        }
        if ((i+1) <= length(nl)) {
            end <- nl[(i+1):length(nl)]
        } else {
            end <- c()
        }
        nl <- c(start, tn[, nl[i]], end)
        idx <- which(! nl %in% labels)
    }

    ord <- sapply(nl, function(i) which(i == labels))

    m <- t(tn)
    nidx <- m %in% labels
    m[nidx] <- -sapply(m[nidx], function(i) which(i== labels))
    m[!nidx] <- sapply(m[!nidx], function(i) unlist(strsplit(i, "_"))[2])
    m <- matrix(as.numeric(m), ncol=2)

    hr <- list(merge=m,
               height=as.double(h),
               order=ord,
               labels=labels,
               method=hclusterResult$method
               )
    stats:::plot.hclust(hr, main=main, xlab=xlab, ylab=ylab, sub="")
}

## set.seed <- 123
## s <- matrix(abs(rnorm(50)), ncol=5)
## rownames(s) <- paste("g", 1:10, sep="_")
## colnames(s) <- paste("t", 1:5, sep="_")
## res <- hcluster(s)
## plotting_hcluster(res)


## perf <- sapply(2:20, function(i) {
##     res <- kmeans(iris[,-5], centers=3, iter.max=i)
##     lev <- sapply(1:3, function(j)
##                   names(which.max(table(iris[res$cluster == j, 5]))))
##     sum(as.numeric(factor(iris[,5], levels=lev)) == res$cluster)/length(iris[,5])

## })
