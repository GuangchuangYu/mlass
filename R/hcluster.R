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

linkage <- function(d, method, cn, cls, dist) {
    sc <- getClosest(d)

    node <- sc$node

    ## if (method == "average") {
    ##     nodes <- node
    ##     ii <- grep("__c", node)
    ##     if (length(ii) > 0) {
    ##         for (i in ii) {
    ##             nodes <- c(nodes, cls[[nodes[i]]]$node)
    ##         }
    ##         nodes <- nodes[-ii]
    ##     }
    ## } else {
    ##     nodes <- node
    ## }

    if ( method == "complete" ) {
        h <- apply(d[node,], 2, max)
    } else if (method == "single") {
        h <- apply(d[node,], 2, min)
    } else {
        ## step 1, expand known clusters
        nodes <- node
        ii <- grep("__c", nodes)
        if (length(ii) > 0) {
            node2 <- nodes[-ii]
        } else {
            node2 <- nodes
        }
        
        while (length(ii) > 0) {
            for (i in ii) {
                node2 <- c(node2, cls[[nodes[i]]]$node)
            }
            ii <- grep("__c", node2)
            nodes <- node2
            if (length(ii) > 0) {
                node2 <- nodes[-ii]
            }
        }
        ## step 2, calcuate mean
        h <- apply(dist[node2,], 2, mean)

        ## step 3, merge known clusters.
        if(length(cls) >=1) {
            for (i in 1:length(cls)) {
                nn <- cls[[i]]$node
                h <- c(h, mean(h[nn]))
                names(h)[length(h)]= names(cls[i])
                h <- h[!names(h) %in% nn]
            }
        }
        h <- h[colnames(d)]
    }

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
    dist <- d
    nr <- nrow(d)
    cls <- list()
    for (i in 1:(nr-1)) {
        cn <- paste("__c", i, sep="-")
        sc <- linkage(d, method, cn, cls, dist)
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
    m[!nidx] <- sapply(m[!nidx], function(i) unlist(strsplit(i, "-"))[2])
    m <- matrix(as.numeric(m), ncol=2)

    hr <- list(merge=m,
               height=as.double(h),
               order=ord,
               labels=labels,
               method=hclusterResult$method
               )
    stats:::plot.hclust(hr, main=main, xlab=xlab, ylab=ylab, sub="")
}



## perf <- sapply(2:20, function(i) {
##     res <- kmeans(iris[,-5], centers=3, iter.max=i)
##     lev <- sapply(1:3, function(j)
##                   names(which.max(table(iris[res$cluster == j, 5]))))
##     sum(as.numeric(factor(iris[,5], levels=lev)) == res$cluster)/length(iris[,5])

## })
