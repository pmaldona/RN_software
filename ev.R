
source("mgg.R")
source("sm.R")

# incremental computation of the evolution potential within a reaction network under additive structural perturbations
# e: current evolution potential (e$rn is the reaction network, e$ev the evolution potential)
# rn: a reaction network (ignored unless e is NULL), s0: starting state (boolean vector)
# M: number of new explorations, n: simulation iterations
ev.evol <- function( e=NULL, rn = if (!is.null(e)) e$rn else rn.testrn(), s0=NULL, M=1, n=5000 ) {
  if (is.null(e)) e <- list(rn=rn,ev=NULL)
  L <- nrow(e$rn$mr)
  if (is.null(s0) && is.null(e$ev)) s0 <- matrix(F,L,1) # starting from the empty set by default
  if (!is.null(s0)) s0[rn.closure(e$rn,s0)] <- T  # the closure of s0
  j <- NULL  # index of the starting state
  if (is.null(e$ev)) {
    e$ev <- list(s=cbind(s0),t=list(list(s0=NULL,sf=NULL))) # s: states, t: their transitions
    rownames(e$ev$s) <- rownames(e$rn$mr)
  }
  else if (!is.null(s0)) {  # the incremental case starting from s0
    j <- which(apply(e$ev$s,2,function(s) all(s==s0)))  # index of s0 in e$ev$s
    if (length(j)==0) {  # s0 is a new state to be added
      e$ev$s <- cbind(e$ev$s,s0)
      e$ev$t <- c(e$ev$t,list(list(s0=NULL,sf=NULL)))
      j <- ncol(e$ev$s)
    }
  }
  if (M>0) for (i in 1:M) {
    repeat {
      if (i>1 || is.null(j))
        j <- sample(ncol(e$ev$s),1)  # the index of the state in e$ev from which the exploration is continued
      s <- e$ev$s[,j]  # the state (logical vector of contained species)
      if (sum(s)<L) break  # s is not complete within rn
    }
    cs <- !s  # complementary species
    k <- sample(sum(cs),1)  # number of added species
    as <- sample(which(cs),k)  # k added species
    s0 <- s  # starting state for the simulation
    s0[as] <- T  # to set added species
    sm <- sm.sim(e$rn,n=n,s0=s0,p=runif(ncol(rn$mr),.3,3),momentum=0,norm=0)
    p <- sm$sa[,n,drop=F]  # last abstract state
    p[rn.closure(e$rn,which(p))] <- T  # its closure
    sf <- which(apply(e$ev$s,2,function(s) all(s==p)))  # index of p in e$ev$s
    if (length(sf)==0) {  # p is a new state to be added
      e$ev$s <- cbind(e$ev$s,p)
      e$ev$t <- c(e$ev$t,list(list(s0=NULL,sf=NULL)))
      sf <- ncol(e$ev$s)
    }
    e$ev$t[[j]]$s0 <- cbind(e$ev$t[[j]]$s0,s0,deparse.level=0)
    e$ev$t[[j]]$sf <- c(e$ev$t[[j]]$sf,sf[1])
  }
  return(e)
}

ev.sort <- function(e) {
  o <- order(colSums(e$ev$s))
  e$ev$s <- e$ev$s[,o]
  e$ev$t <- e$ev$t[o]
  for (i in 1:length(e$ev$t))
    if (!is.null(e$ev$t[[i]]$sf))
      e$ev$t[[i]]$sf <- o[e$ev$t[[i]]$sf]
  mt <- matrix(0,ncol(e$ev$s),ncol(e$ev$s))  # transition matrix
  for (i in 1:length(e$ev$t))
    if (!is.null(e$ev$t[[i]]$sf))
      for (j in e$ev$t[[i]]$sf) mt[j,i] <- mt[j,i] + 1
  e$mt <- mt
  return(e)
}
