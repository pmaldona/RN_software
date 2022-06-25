
# state and flow perturbations

# randomize components of a state or flow vector using a log normal distribution
# v vector, mask of affected components (defaults to active ones), mu and sigma parameters of the distribution
pert.randomize <- function(v,mask=(v!=0),mu=0,sigma=.5) {
  mask <- rep(T,length(v)) & mask
  k <- which(mask)
  if (length(k)==0) return(v)
  v[k] <- exp(rnorm(length(k),mu,sigma))
  return(v)
}

# activation of components of v (changing values from 0 to 1 or from >0 to 0), mask (defaults to all components) 
# n total active components required, up to p active components to be preserved (defaults to all active components)
pert.activation <- function(v,mask=rep(T,length(v)),n=sum(mask),p=sum(v[mask]>0)) {
  w <- v[mask]
  l <- length(w)
  a <- sum(w>0)
  n <- max(min(n,l),0)
  p <- max(min(p,a),0)
  if (n==0) w <- 0*w
  else if (n<=p) {
    if (a==1) k <- which(w>0) else k <- sample(which(w>0),n)
    w[-k] <- 0
  }
  else if (p==0) {
    if (l==1) k <- 1 else k <- sample(l,n)
    w[-k] <- 0
    w[k[w[k]==0]] <- 1
  }
  else {  # n>p and p>0
    i <- which(w>0)
    if (length(i)>1) i <- sample(i,p)  # the ones we keep
    j <- (1:l)[-i]  # the remaining components
    k <- sample(length(j),n-p)
    w[j[-k]] <- 0
    k <- j[k]
    w[k[w[k]==0]] <- 1
  }
  v[mask] <- w
  return(v)
}

# perturb a vector by adding or substracting components and randomizing the result
# v vector, d component delta, at least nmin active components, sigma random dispersion
pert.delta <- function(v,d=1,nmin=5,sigma=.5) {
  n <- max(sum(v>0)+d,nmin)
  v <- pert.activation(v,n=n)
  v <- pert.randomize(v,sigma=sigma)
  return(v)
}

# creates a structure to study evolutive paths of a reaction network rn
# result: the reaction network rn, a random walk list rw
# each rw[[i]] contains matrices f, s, p, c, a and r were each column stores
# the starting state (f,s), the perturbed state p, its convergence (c), its abstraction (a) and used species (u)
# for each perturbation and simulation step applied
pert.start <- function(rn=rn.testrn()) {
  s <- rep(0,nrow(rn$mr))
  return(list( rn=rn,species=rownames(rn$mr),reactions=colnames(rn$mr),
               rw=list(list(f=NULL,s=NULL,p=NULL,c=NULL,a=NULL,u=NULL ))))
}

# simulates the dynamics from state s with flow vector f and returns the end result
# cutoff is the concentration threshold for species to be reactive
# n is the maximum number of iterations (the simulation can stop before)
pert.simul <- function(rn,s,f,cutoff=.1,n=5000) {
  sm <- sm.maksim(rn,n=n,s0=s,p=f,e=cutoff,inflow=F)
  if (!is.null(ncol(sm$s))) return(sm$s[,ncol(sm$s)])
  else return(NULL)
}

# returns the boolean abstract state (closure) for a given reaction network rn
# species concentrations vector s and flow constants vector f
pert.abstract <- function(rn,s,f,cutoff=.1) {
  k <- rn.closure(rn,which(s>=cutoff),which(f>0))
  v <- rep(F,length(s))
  v[k] <- T
  return(v)
}

# returns the boolean vector of used species for a given reaction network rn
# species concentrations vector s and flow constants vector f
pert.used <- function(rn,s,f,cutoff=.1) {
  k <- rn.used(rn,which(s>=cutoff),which(f>0))
  v <- rep(F,length(s))
  v[k] <- T
  return(v)
}
