#' Calculate the probabilities of occurrence for every species at every iteration step.
#' 
#' Here the probability of occurrence is considered simply as the fraction of iteration steps
#' where a species is present at the end of the simulation
#' 
#' @param e
#'
#' @return p, a matrix such that
#' p_i_j is the number of times that species i was present at the end of the simulation of the first j iteration steps
occurrence_probabilities <- function(m){
  c <- matrix(nrow=nrow(m) , ncol=ncol(m))
  for (i in 1:nrow(c)){
    c[i,1] <- m[i,1]
    for (j in 2:ncol(c)){
      c[i,j]<-c[i,j-1]
      if(m[i,j] == 1){
        c[i,j]<-c[i,j]+1
      }
    }
  }
  for (i in 1:nrow(c)){
    for (j in 2:ncol(c)){
      c[i,j]<-c[i,j]/j
    }
  }
  levelplot(c)
  return(c)
}