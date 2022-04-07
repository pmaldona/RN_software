library(lattice)

#' A method to calculate the conditional probabilities of species occurrences
#'
#' @param abstraction
#' A numeric matrix where each column represents an abstraction of a state of a dynamical system
#' @param m 
#' The index of the species thats presence is the condition we are interested in
#'
#' @return P_c
#' A matrix such that P_c[i,j] is the probability that species i is present at timestep j given
#' that m is also present
conditional_probabilities <- function(abstraction, m){
  M_0 <- matrix (0, nrow(abstraction), ncol(abstraction))
  M_1 <- matrix (0, nrow(abstraction), ncol(abstraction))
  P_c <- matrix (0, nrow(abstraction), ncol(abstraction))
  for(j in 1:ncol(abstraction)){
    for(i in 1:nrow(abstraction)){
      if(j>1){
      M_0[i,j] <- M_0[i,j-1]
      M_1[i,j] <- M_1[i,j-1]
      }
      if(abstraction[m,j]==1){
        if(abstraction[i,j]==0){
          M_0[i,j] <- M_0[i,j]+1
        }
        else{
          M_1[i,j] <- M_1[i,j]+1
        }
      }
      if((M_0[i,j]+M_1[i,j])>0){
        P_c[i,j] <- M_1[i,j] / (M_0[i,j]+M_1[i,j])
      }
    }
  }
  return(P_c)
}
