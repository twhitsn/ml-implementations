nb.classify <- function(form, df){
    dm <- model.frame(form, df)
    
    classes <- dm[,1]
    train <- dm[,-1]

    matrices <- list()
    priors <- list()
    
    for(uc in unique(classes)){
        priors[[uc]] <- list()
        
        for(ci in 1:ncol(train)){
            column.name <- colnames(train)[ci]
            priors[[uc]][[column.name]] <- list()
            column <- train[,ci]
            
            # gaussian, will predict later with dnorm
            if(is.numeric(column)){
                priors[[uc]][[column.name]]['mean'] <- mean(column)
                priors[[uc]][[column.name]]['sd'] <- sd(column)
                
            # count
            } else {
                uv <- unique(column)
                
                for(ri in 1:length(uv)){
                    cur.value <- uv[ri]
                    priors[[uc]][[column.name]][[toString(cur.value)]] <- length(which(column == cur.value & classes == uc)) / length(column)
                }
            }
        }
    }
    
    return(priors)
}

nb.predict <- function(fit, test){
    pred <- rep(0, nrow(test))
    for(i in 1:nrow(test)){
        row <- test[i,]
        
        probs <- sapply(1:length(row), function(j){
            name <- names(row)[j]
            value <- row[[j]]
            
            sapply(fit, function(l){
                if(is.numeric(value)){
                    dnorm(value, l[[name]][['mean']], l[[name]][['sd']])
                } else {
                    
                    # need toString, or value is considered a factor
                    l[[name]][[toString(value)]]
                }
            })
        })
        
        prods <- apply(probs, 1, prod)
        pred[i] <- names(prods)[which.max(prods)]
    }
    
    return(pred)
}