###############################################
#
#          Written by Tim Whitson
#
###############################################

kmeans <- function(data, K, pp = FALSE, plot.data = FALSE, iter.max = 10){
    data <- data.matrix(data)
    
    if(K > nrow(data)){
        stop("Error: Cannot have more clusters than data points.")   
    }
            
    # data plotting
    plot_data <- function(centroids, labels){
        plot.new()
        
        # initial plot
        if(length(labels) == 0){
            plot(data)
            
        # final plot
        } else {
            plot(data)
            # create random list of colors, length K
            colors = sample(colours(), K)
            
            # plot each cluster
            for(i in 1:K){
                cluster <- get_cluster(i, labels)
                points(cluster, col = colors[i], pch = 19)
            }
            
            # plot centroids
            points(centroids, col = colors, pch = 3)
        }
    }
        
    # get euclidean distance of each point from centroids
    get_distances <- function(centroids){
        distances <- sapply(1:K, function(k){
            sapply(1:nrow(data), function (i){
                dist(rbind(data[i,], centroids[k,]))
            })
        })
        
        return(distances)
    }
        
    # match each point to closest centroid
    get_labels <- function(distances){
        new_labels <- apply(distances, 1, function(d){
            which.min(d)
        })
        return(new_labels)
    }
        
    # calculate new centroids by mean of data points
    update_centroids <- function(labels){
        return(matrix(
            sapply(1:K, function(k) colMeans(get_cluster(k, labels))), 
            ncol = ncol(data),
            byrow = TRUE
        ))
    }
        
    # get cluster of data based on indexes of cluster labels
    # drop == FALSE, otherwise if only one row, colMeans will return error
    get_cluster <- function(index, labels){
        cluster <- data[which(labels == index),, drop = FALSE]
        return(cluster)
    }
        
    # check for convergence
    get_sse <- function(distances){
        new_sse <- sum(apply(distances, 1, function(d){
            min(d) ^ 2
        }))

        return(new_sse)
    }

    # iterate, recomputing centroids and finding nearest centroids for each point
    centroids <- data[sample(1:nrow(data), K),]
    distances <- get_distances(centroids)
    labels <- get_labels(distances)
    sse <- get_sse(distances)
    
    iterations <- 0

    while(iterations < iter.max){
        iterations <- iterations + 1

        # maintain K clusters
        while(length(unique(labels)) != K){
            missing <- head(setdiff(1:K, unique(labels)), 1)
            farthest <- data[which.max(apply(distances, 1, function(d){
                max(d)
            })),]

            centroids[missing,] <- farthest
            
            distances <- get_distances(centroids)
            labels <- get_labels(distances)
        }
        
        centroids <- update_centroids(labels)
        distances <- get_distances(centroids)
        new_labels <- get_labels(distances)
        
        sse <- get_sse(distances)
        
        if(all(labels == new_labels)){
            break   
        } else {
            labels <- new_labels   
        }
    }

    if(iterations >= iter.max){
        warning(paste('Did not converge within maximum number of iterations =', iter.max))
    }
    
    if(plot.data == TRUE){
        plot_data(centroids, labels)   
    }
    
    return(list(
        'centroids' = centroids,
        'labels' = labels,
        'iterations' = iterations,
        'sse' = sse
    ))
}