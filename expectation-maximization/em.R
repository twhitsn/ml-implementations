###################################################################
#
#  Expectiation-Maximization Algorithm for Gaussian Distribution
#
#                         Tim Whitson
#
###################################################################

library(mvtnorm)

em <- function(data, K, iter.max = 100, epsilon = 1e-5){
    data <- data.matrix(data)
    dimensions <- dim(data)
    
    # initialize parameters - 
    #   mu      = K random points (or mu input)
    #   sigma   = identity matrix, d x d
    #   prior   = 1 / K
    init <- function(){
        params = list()
        
        params$u <- data[sample(1:dimensions[1], K),]
        
        params$s <- lapply(1:K, function(i){
            diag(dimensions[2])
        })
        
        params$p <- sapply(1:K, function(i){
            1 / K  
        })
        
        params$l <- log.likelihood(params)
        
        return(params)
    }
    
    # log likelihood
    log.likelihood <- function(params){
        l <- sum(log(rowSums(sapply(1:K, function(i){
            params$p[i] * ns.dmvnorm(params$u[i,], params$s[[i]])
        }))))
        
        return(l)
    }
    
    # e-step
    expectation <- function(params){
        w <- c()
        
        # probability density
        for(i in 1:K){
            w <- cbind(w, params$p[i] * ns.dmvnorm(params$u[i,], params$s[[i]]))
        }
        
        # divide each value by its row sum
        w <- sweep(w, 1, rowSums(w), '/')
        
        return(w)
    }
    
    # m-step
    maximization <- function(w){
        new_params <- list()
        
        u <- matrix(nrow = K, ncol = dimensions[2])
        s <- c()
        
        for(i in 1:K){
            w.i <- w[,i]
            
            # update means
            u[i,] <- colSums(data * w.i) / sum(w.i)
            
            # update variance
            dif <- sweep(data, 2, u[i,])
            new.s <- t(dif) %*% (dif * w.i) / sum(w.i)
            s[[i]] <- matrix(new.s, ncol = dimensions[2], byrow = TRUE)
        }
        
        new_params$u <- u
        new_params$s <- s
        
        # update probability
        new_params$p <- colSums(w) / dimensions[1]
        
        return(new_params)
    }
    
    # dmvnorm singular matrix fix
    # add .001 to diagonal
    ns.dmvnorm <- function(mu, sigma){
        dn <- dmvnorm(data, mu, sigma)

        if(all(dn == 0)){
            new.s <- sigma
            diag(new.s) <- diag(new.s) + .001
            dn <- dmvnorm(data, mu, new.s)
        }
        
        return(dn)
    }
    
    # run EM
    run <- function(){
        params <- init()
        iter.n <- 0
        
        while(iter.n < iter.max){
            # e-step
            w <- expectation(params)
            
            # iterate here, so if cluster is lost and has to be reset, does not count as iteration
            iter.n <- iter.n + 1
            
            # m-step
            new_params <- maximization(w)
            
            # distance between means
            u.dist <- sum(sqrt(rowSums((params$u - new_params$u) ^ 2)))
            
            # log likelihood
            new_params$l <- log.likelihood(new_params)
            
            converge <- (new_params$l - params$l) < epsilon
            params <- new_params
            
            if(converge){
                break
            }
        }
        
        labels <- apply(w, 1, function(c){
            which.max(c)
        })
        
        # return values
        params$w <- w
        params$labels <- labels
        params$iterations <- iter.n
        
        return(params)
    }
    
    return(run())
}