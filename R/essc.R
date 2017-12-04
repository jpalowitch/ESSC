#' essc
#'
#' Identify statistically significant communities in undirected networks
#' @param Adj.Matrix Adjacency matrix of the network for which you'd like to find communities
#' @param alpha False discovery rate. Value between 0 and 1 that sets an upper limit on the false discovery rate of the hypothesis testing procedure in ESSC
#' @param Null The null distribution used for comparing the observed number of connections between a single vertex and a community of vertices. Can be set to either "Binomial" or "Poisson." Default is "Binomial"
#' @param Num.Samples The number of randomly selected neighborhoods used to initiate the ESSC algorithm. The maximum number is the number of vertices in Adj.Matrix. Defaults to the number of vertices in Adj.Matrix
#' @keywords community detection, community extraction
#' @return a list containing the objects
#' \itemize{
#'    \item Communities: a list of identified communities
#'    \item Background: a numeric vector identifying vertices that did not belong to any statistically significant community
#' }
#'@references
#'\itemize{
#'     \item Wilson, James D.,Wang, Simi, Mucha, Peter J., Bhamidi, Shankar, and Nobel, Andrew B. (2014). “A testing based extraction algorithm for identifying significant communities in networks.”
#'     The Annals of Applied Statistics Vol. 8, No. 3, 1853-1891
#' }
#' @author James D. Wilson
#' @examples
#' net <- stochastic.block(n = 1000, k = 3, P = cbind(c(0.1, 0.01, 0.01), c(0.01, 0.1, 0.01), c(0.01, 0.01, 0.1)), sizes = c(300, 300, 400))
#' results <- essc(net$Adjacency, alpha = 0.10, Null = "Poisson")
#' print(results$Communities)
#' @export

###essc Function. Run this function to extract communities
essc = function(Adj.Matrix, alpha, Null = c("Binomial", "Poisson"),
                Num.Samples = nrow(Adj.Matrix), max.iter=300, verbose=FALSE,
                tau=0.9, exhaustive=TRUE){
    symdiff = function(X,Y){setdiff(union(X,Y),intersect(X,Y))}

    # Function to filter overlap
    filter_overlap <- function (comms, tau) {

      K <- length(comms)
      scores <- unlist(lapply(comms, length))

      jaccard_mat0 <- matrix(0, K, K)
      for (i in 1:K) {
        for (j in 1:K) {
          jaccard_mat0[i, j] <- length(intersect(comms[[i]], comms[[j]])) /
            length(comms[[i]])
        }
      }

      jaccard_mat <- jaccard_mat0
      diag(jaccard_mat) <- 0
      max_jacc <- max(jaccard_mat)
      deleted_comms <- integer(0)

      while (max_jacc > tau) {

        inds <- which(jaccard_mat == max_jacc, arr.ind = TRUE)[1, ]

        # keep smaller comm
        delete_comm <- inds[which.min(c(scores[inds[1]], scores[inds[2]]))]
        jaccard_mat[delete_comm, ] <- 0
        jaccard_mat[, delete_comm] <- 0
        deleted_comms <- c(deleted_comms, delete_comm)
        max_jacc <- max(jaccard_mat)

      }

      kept_comms <- setdiff(1:K, deleted_comms)

      return(list("final_comms" = comms[kept_comms],
                  "kept_comms" = kept_comms))

    }

    # Set-up
    degrees <- rowSums(Adj.Matrix)
    n <- length(degrees)
    index <- 1:dim(Adj.Matrix)[1]
    extractFrom <- sample(index,Num.Samples,replace = FALSE)

    #initialize
    Community <- rep(list(NULL),Num.Samples)
    unq <- TRUE
    nodes <- matrix(0,2,Num.Samples)
    nodes[1, ] <- extractFrom
    which_match <- 0


    i <- 1
    count <- 0
    clustered_nodes <- integer(0)
    #Running long loop
    for(node in extractFrom){

      if (node %in% clustered_nodes)
        next

        count <- count+1
        #cat("######################################################\n")
        #cat(paste0("extraction ",count,"\n"))
        B0 <- c(node,which(Adj.Matrix[node,]>0))
        if (verbose) {
          cat("###### Starting Main.Search Loop number",
              match(node, extractFrom), "######\n")
        }
        temp <- Main.Search(Adj.Matrix, alpha, B0, Null, max.iter, verbose)

        if(length(temp$Community) > 1){
            if(i > 1){
                unq <- TRUE
                which_match <- 0
                for(j in 1: i){
                    if(length(symdiff(Community[[j]],temp$Community)) == 0){
                        unq <- FALSE
                        which_match <- j
                    }
                }
            }
            if(unq){
                Community[[i]] <- temp$Community
                nodes[2, count] <- i
                i = i+1
            } else{
                nodes[2, count] <- which_match
            }
        } else{
            nodes[2, count] <- 0
        }

      clustered_nodes <- unlist(Community)
    }

    Community <- Community[unlist(lapply(Community,function(comm)length(comm)>0))]

    if(length(Community) == 0)
        Community <- list(NULL)

    if (length(Community) > 1) {
      Community <- filter_overlap(Community, tau)$final_comms
    }

    Background <- setdiff(1:n,unlist(Community))

    return(list(Communities = Community, Background = Background))
}

###Single Search Function

Main.Search = function(Adj.Matrix,alpha,B0, Null, max.iter=300, verbose=FALSE){
    j <- 0
    degrees <- rowSums(Adj.Matrix)
    n <- length(degrees)
    #Ensuring B1 and B0 are not equal to start the loop
    B1 <- integer(0)
    #Function for checking equality of matrices
    vectorequal = function(x,y){
        is.numeric(x) && is.numeric(y) && length(x) == length(y) && all(x==y)
    }
    while(vectorequal(B0,B1) == FALSE & j <= 30){

        j <- j + 1
        if(verbose && (j < 5 || j%%5==0)) {
           cat(paste("iteration",j,"\n"))
        }
        if(j > 1)
            B0 <- B1
	  if(length(B0)>1)
	        duBs <- rowSums(Adj.Matrix[, B0])
	  if(length(B0)==1)
		  duBs <- Adj.Matrix[, B0]
        pB <- sum(degrees[B0])/sum(degrees) #probability of connection to B
        if(Null == "Binomial")
            pvals <- pbinom(duBs,degrees,pB,lower.tail = FALSE)
        if(Null == "Poisson")
            pvals <- ppois(duBs,degrees*pB,lower.tail = FALSE)
        pvals[degrees==0] <- 1
        pvals_bh <- pvals*n/rank(pvals)
        if(sum(pvals_bh <= alpha) == 0)
        {
          B1 <- integer(0)
          if (verbose) {
            cat("B1 reached size 0 in Main.Search\n")
          }
          break
        }
        threshold <- max(pvals[pvals_bh<=alpha])
        B1 <- which(pvals<=threshold)
        if(length(B1) <1 ){break}
    }
    Community <- B1
    if(j > max.iter){
        if (verbose) {
          cat("Main.Search failed to converge before 30 iterations\n")
        }
        Community <- integer(0)
    }
    return(list(Community = Community))
}


