
source("mgg.R")
source("sm.R")
source("pert.R")

# distribution functions for the size of generated perturbations (n is the maximum perturbation size)
ev.Pexp <- function(n) exp(-(1:n))
ev.P <- function(n) (n+1)/2 - abs((n+1)/2 - (1:n))

# incremental computation of the evolution potential within a reaction network under additive structural perturbations
# e: current evolution potential (e$rn is the reaction network, e$ev the evolution potential)
# rn: a reaction network (ignored unless e is NULL), s0: starting state (boolean vector), P: perturbation size distribution
# M: number of new explorations, n: simulation iterations, systematic: always s0 as starting point vs selecting randomly
ev.evol <- function( e=NULL, rn = if (!is.null(e)) e$rn else rn.testrn(), s0=NULL, P=ev.P, M=1, n=5000, systematic=T ) {
  if (is.null(e)) e <- list(rn=rn,ev=NULL)
  rSums <- colSums(rn$mr) # reactants sums for each reaction (used to determine active reactions given a state)
  L <- nrow(e$rn$mr) # number of species in the reaction network
  if (is.null(s0) && (is.null(e$ev)|systematic)) s0 <- matrix(F,L,1) # starting from the empty set by default
  if (!is.null(s0)) s0[rn.closure(e$rn,s0)] <- T  # the closure of s0
  j <- NULL  # index of the starting state
  if (is.null(e$ev)) {
    e$ev <- list(s=cbind(s0),a=0,t=list(list(s0=NULL,sf=NULL))) # s: states, a: attractiveness, t: their transitions
    rownames(e$ev$s) <- rownames(e$rn$mr)
  }
  if (!is.null(s0)) {  # the incremental case starting from s0
    j <- which(apply(e$ev$s,2,function(s) all(s==s0)))  # index of s0 in e$ev$s
    if (length(j)==0) {  # s0 is a new state to be added
      e$ev$s <- cbind(e$ev$s,s0)
      e$ev$a <- c(e$ev$a,0)
      e$ev$t <- c(e$ev$t,list(list(s0=NULL,sf=NULL)))
      j <- ncol(e$ev$s)
    }
  }
  if (M>0) for (i in 1:M) {
    repeat {
      if (!systematic & (i>1 || is.null(j)))
        j <- sample(ncol(e$ev$s),1)  # the index of the state in e$ev from which the exploration is continued
      s <- e$ev$s[,j]  # the state (logical vector of contained species)
      if (sum(s)<L) break  # ok if s is not the whole network (if it is, we choose another starting state)
    }
    rS <- colSums(e$rn$mr[s,,drop=F])  # reactants sums
    nr <- rS!=rSums  # non active reactions
    k <- sample(sum(nr),1,prob=P(sum(nr)))  # number of added reactions
    ar <- sample(which(nr),k)  # k added reactions
    as <- rowSums(e$rn$mr[,ar,drop=F])>0  # added species (reactants of added reactions)
    s0 <- s  # starting state for the simulation
    s0[as] <- T  # to set added species
    s0[rn.closure(e$rn,which(s0))] <- T  # the closure of s0
    sm <- sm.maksim(e$rn,n=n,s0=s0,p=runif(ncol(rn$mr),.3,3))
    # sm <- sm.sim(e$rn,n=n,s0=s0,p=runif(ncol(rn$mr),.3,3),e=c(.1,.1),w=Inf,momentum=0)
    # sm <- sm.sim(e$rn,n=n,s0=s0,p=runif(ncol(rn$mr),.3,3),momentum=0)
    # sm <- sm.cfsim(e$rn,n=n,s0=s0)
    p <- sm$sa[,ncol(sm$sa),drop=F]  # last abstract state
    p[rn.closure(e$rn,which(p))] <- T  # its closure
    sf <- which(apply(e$ev$s,2,function(s) all(s==p)))  # index of p in e$ev$s
    if (length(sf)==0) {  # p is a new state to be added
      e$ev$s <- cbind(e$ev$s,p)
      e$ev$a <- c(e$ev$a,0)
      e$ev$t <- c(e$ev$t,list(list(s0=NULL,sf=NULL)))
      sf <- ncol(e$ev$s)
    }
    e$ev$t[[j]]$s0 <- cbind(e$ev$t[[j]]$s0,s0,deparse.level=0)
    e$ev$a[sf[1]] <- e$ev$a[sf[1]] + 1
    e$ev$t[[j]]$sf <- c(e$ev$t[[j]]$sf,sf[1])
  }
  return(e) # e$rn reaction network, e$ev$s end (or starting) states, e$ev$t[[i]] transitions generated from e$ev$s[,i]
            # e$ev$t[[i]]$s0 perturbed states generated, e$ev$t[[i]]$sf final states reached simulating the dynamics
}

# analysis of the results of the evolution
ev.analyse <- function(e) {
  o <- do.call("order",unname(as.data.frame(t(e$ev$s))))
  o <- o[order(colSums(e$ev$s)[o])]
  e$ev$s <- e$ev$s[,o]
  e$ev$a <- e$ev$a[o]
  oi <- o; oi[o] <- 1:length(o)
  for (i in 1:length(e$ev$t)) e$ev$t[[i]]$sf <- oi[e$ev$t[[i]]$sf]
  s0 <- NULL; sf <- NULL  # initial (s0, boolean vector of species) and final (sf, index to s) states
  for (i in 1:length(e$ev$t)) { s0 <- cbind(s0,e$ev$t[[i]]$s0); sf <- c(sf,e$ev$t[[i]]$sf) }
  o <- do.call("order",unname(as.data.frame(t(s0))))
  o <- o[order(colSums(s0)[o])]
  s0 <- s0[,o]; sf <- sf[o]
  p <- matrix(0,nrow(s0),0)  # perturbed states (empty at the start, unique s0's at the end)
  tpf <- matrix(0,ncol(e$ev$s),0)  # transition matrix from perturbed states to final states (empty at the start)
  k <- 0
  for (i in 1:ncol(s0))
    if (i==1 || !all(s0[,i]==s0[,i-1])) {  # a new perturbed state (different from the preceding one)
      k <- k + 1; p <- cbind(p,s0[,i]); tpf <- cbind(tpf,0); tpf[sf[i],k] <- 1
    }
    else {
      tpf[sf[i],k] <- tpf[sf[i],k] + 1
    }
  e$p <- p; e$tpf <- tpf
  dpf <- sapply(1:ncol(p), function(i) sapply(1:ncol(e$ev$s), function(j) sum(abs(p[,i]-e$ev$s[,j]))))
  e$dpf <- dpf  # distance from perturbed states to final states
  cs <- colSums(tpf)
  ntpf <- sweep(tpf,2,cs+(cs==0),"/")  # column normalized transition matrix
  pdf <- t(sapply(0:nrow(e$rn$mr), function(i) rowSums(dpf==i)))  # perturbations at given distances from final states
  cpdf <- t(sapply(0:nrow(e$rn$mr), function(i) rowSums((dpf==i)*ntpf)))  # the ones that converge to the final states
  e$ntpf <- ntpf; e$pdf <- pdf; e$cpdf <- cpdf
  m <- cpdf/(pdf+(pdf==0))
  e$lat <- colSums(sweep(m,1,1:nrow(pdf),"*")*(pdf>0))
  e$prec <- colSums(sweep(1-m,1,1/(1:nrow(pdf)),"*")*(pdf>0))
  e$atr <- colSums(cpdf); e$atr <- e$atr/sum(e$atr)
  return(e)
}

ev.disp <- function(ae,z=1,save=F) {
  m <- ae$cpdf/ae$pdf
  g <- mgg.img(m,ylim=c(0,nrow(ae$p)),zlim=c(0,z),col=viridis(256)) +
    labs( title="Probability of converging to a state from a distance",
          x="state (ordered by number of species)",y="distance" )
  print(g); if (save) mgg.save("prob.png",dim=c(7,3.5))
  g <- mgg.lin(ae$atr) +
    labs(x="state (ordered by number of species)",y="attractiveness (probability of convergence)") +
    theme(legend.position="none") 
  print(g); if (save) mgg.save("atract.png",dim=c(7,3.5))
  g <- mgg.lin(rbind(latitude=ae$lat,precariousness=ae$prec,resistance=ae$lat/ae$prec,by_size=sort(ae$lat/ae$prec,decreasing=T))) +
    labs(x="state (ordered by number of species)",y=NULL)
  print(g); if (save) mgg.save("lat_prec.png",dim=c(7,3.5))
  g <- mgg.scat(rbind(latitude=ae$lat,precariousness=ae$prec)) +
    labs(x="latitude",y="precariousness") + theme(legend.position="none")
  print(g); if (save) mgg.save("prec_vs_lat.png",dim=c(7,3.5))
}

ev.do <- function(sm=T) {
  mP <- rbind(c(0,0,0,1,1,50),c(0,0,0,1,2,20),c(0,0,0,2,1,20),c(0,0,0,2,2,10))
  Nr <- 100
  S <- paste0("s",1:Nr)
  C <- 1:3
  P <- function() rg.rsrt(m=mP)
  D <- list(s=length(S)/(1:length(S))+length(S)/5,c=length(C)/(1:length(C)))
  rn <<- rg.rn(Nr=Nr,S=S,C=C,P=P,D=D)
  trn <<- rn.trim(rn)
  if (sm) {
    e <<- ev.evol(rn=trn,M=5000,n=1000)
    ae <<- ev.analyse(e)
    ev.disp(ae)
  }
}
