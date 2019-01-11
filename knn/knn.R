kNN <- function(form, train.df, test, k, dist.method = 'euclidean'){
    train.dm <- model.frame(form, train.df)
    
    train.class <- train.dm[,1]
    train <- train.dm[,-1]
    
    mode <- function(v){
        names(tail(sort(table(v)), 1))
    }
    
    # get all distances from test points
    distances <- sapply(1:nrow(test), function(i){
        sapply(1:nrow(train), function(j){
            dist(rbind(as.numeric(test[i,]), as.numeric(train[j,])), method = dist.method)
        })
    })
    
    # get k closest distances (by index)
    min.distances <- apply(distances, 2, function(i){
        order(i)[1:k]
    })
    
    majority <- apply(min.distances, 2, function(i){
        mode(train.class[i])        
    })
    
    return(majority)
}
