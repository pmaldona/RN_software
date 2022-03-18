
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
# result: the reaction network rn, current reaction constants and species concentrations f and s,
# hf, hs and hc are matrices were each column stores the starting state (f,s) and the concentration convergence (c)
# for each dynamics simulation applied
pert.start <- function(rn=rn.testrn()) {
  s <- rep(0,nrow(rn$mr))
  f <- pert.randomize(rep(1,ncol(rn$mr)))
  return(list(rn=rn,f=f,s=s,hf=NULL,hs=NULL,hc=NULL))
}

# the variable v, "s" (species concentrations, default) or "f" (reaction constants), is perturbed by function op
# the current state (f,s) is modified according to the perturbation applied
# several perturbation can be applied sequentially to modify the current state
pert.apply <- function(e,v="s",op=pert.delta,...){
  e[[v]] <- op(e[[v]],...)
  return(e)
}

# simulates the dynamics from the current state (after applying perturbations) and stores the result in h matrices
pert.simul <- function(e,cutoff=.1,n=5000) {
  sm <- sm.maksim(e$rn,n=n,s0=e$s,p=e$f,e=cutoff,inflow=F)
  e$hf <- cbind(e$hf,e$f)
  e$hs <- cbind(e$hs,e$s)
  if (!is.null(ncol(sm$s))) e$s <- sm$s[,ncol(sm$s)]
  e$s[e$s<cutoff] <- 0
  e$hc <- cbind(e$hc,e$s)
  return(e)
}

# returns the boolean matrix of abstract states (closure) reached by simulating the dynamics
pert.abstract <- function(e,cutoff=.1) {
  if (is.null(e$hc)) return(NULL) # no history yet, NULL result
  closure <- function(i) {
    s <- rn.closure(e$rn,which(e$hc[,i]>=cutoff),which(e$hf[,i]>0))
    v <- rep(F,nrow(e$hc))
    v[s] <- T
    return(v)
  }
  return(sapply(1:ncol(e$hc),closure))
}
