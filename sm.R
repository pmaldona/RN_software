
source("rn.R")
source("rg.R")

# generates random networks of L reactions according to s model
sm.genrn <- function(L=12,s=1) {
  lr <<- switch(s,
                rg.lreac(L=L,U=L,g=rg.genreac(dCM=rg.dCM(C=1:2,M=1))), # (1) unbalanced
                rg.lreac(L=L,g=rg.genreac(frt=function() rg.rrt(pmix=.5))), # (2) with balanced mass and lots of autocatalytic react.
                rg.lreac(L=L,U=L,g=rg.genreac(frt=function() rg.rrt(pmix=.5))), # (3) the same unbalanced
                rg.lreac(L=L,U=L,g=rg.genreac(dCM=rg.dCM(M=1),frt=function() rg.rrt(N=c(3,3)))), # (4) unbalanced with up to 3 reac./prd.
                rg.lreac(L=L,U=L,g=rg.genreac(dCM=rg.dCM(C=1:2,M=1,dC=c(1000,1)),frt=rg.rsrt)), # (5) unbalanced
                rg.lreac(L=L,g=rg.genreac(dCM=rg.dCM(C=1:2,M=1,dC=c(1000,1)),frt=rg.rsrt)), # (6) balanced
                rg.lreac(L=L,g=rg.genreac(dCM=rg.dCM(C=1,M=1), frt=function() rg.rsrt(m=cbind(0,0,0,1,1,1)))) # (7) balanced
  )
  ucn0 <<- rg.ucnet(lreac=lr)
  ucn <<- ucn0$copy()
  rg.cyclcon(ucn)
  rg.rcon(ucn,nrow(ucn$mp)/2)
  rn <<- ucn$rn()
}

# Euler method simulation of a reaction network rn for n iterations at dt time step (adaptative)
# with initial state s0 (infused up to time t0), mass action parameter p (vector, one component per reaction),
# minimum activation threshold e (vector: low, high), maximum species concentration w, acceptance of deficit of species,
# momentum (0 => no momentum), and normalization norm (norm==0 => no normalization).
sm.sim <- function( rn,n=1000,dt=.1,s0=runif(nrow(rn$mr),.9,1.1),t0=dt*n/10,p=runif(ncol(rn$mr),.9,1.1),
                    e=1e-3*c(1,10), w=1e1, deficit=T, momentum=.95, norm=0, inflow=T ) {
  s <- matrix(0,nrow(rn$mr),n) # state matrix (species concentration, one column vector per iteration)
  sa <- matrix(F,nrow(rn$mr),n) # state matrix abstraction
  rownames(s) <- rownames(rn$mr)
  v <- matrix(0,ncol(rn$mr),n) # process matrix (reactions activity, one column vector per iteration)
  rownames(v) <- colnames(rn$mr)
  vm <- numeric(ncol(rn$mr)) # process vector with momentum
  t <- numeric(n) # time vector
  m <- rn$mp - rn$mr
  if (length(e)<2) e <- c(e,e)
  cs <- logical(nrow(rn$mr)) # current state (species available)
  if (inflow==F)
    inflow <- numeric(0)
  else if (inflow==T)
    inflow <- which(rowSums(rn$mp[,which(colSums(rn$mr)==0),drop=F])>0)
  # print(inflow)
  f <- 1 # time dilation factor
  for (i in 1:n) {
    if (i>1) {
      ds <- m %*% cbind(vm) # change of concentrations per time unit
      k <- which(ds<0) # species that are losing concentration
      f <- 1 # factor to be applied to ds
      if (norm>0 && length(k)>0) { # normalization (only species that are losing concentration)
        f <- 1/max(-ds[k])
      }
      if (!deficit)
        f <- min(-s[k,i-1]/ds[k]/dt,f) # f may need to be shorter to avoid deficit
      t[i] <- t[i-1] + f*dt
      s[,i] <- pmin(s[,i-1]+f*ds*dt,w)
    }
    if (t[i]<=t0) s[,i] <- s[,i] + s0*(if (t0<=0) 1 else f*dt/t0)
    if (length(inflow)>0) s[inflow,i] <- s0[inflow]
    cs[s[,i]>e[2]] <- T # species over high threshold are available
    cs[s[,i]<=e[1]] <- F # species below low threshold are not available
    sa[,i] <- cs
    k.s <- which(cs) # available species
    k.r <- which(rbind(cs==F) %*% rn$mr == 0) # reactions with all reactants available
    if (length(k.r)==0) next
    v[k.r,i] <- dt*p[k.r]*exp(rbind(log(s[k.s,i])) %*% rn$mr[k.s,k.r])
    vm <- momentum*vm + (1-momentum)*v[,i]
  }
  return(list(s=s,sa=sa,v=v,t=t))  # s: species concentration, sa: species abstract existence, v: flow vector, t: time
}

# returns the sequence of state abstractions from de simulation results sm of reaction network rn
sm.abstr <- function(rn,sm) {
  l <- ncol(sm$v)
  ia <- c(1, 1 + which(colSums(sm$sa[,2:l]!=sm$sa[,1:(l-1)])>0)) # starting index of available species state bouts
  v <- sm$sa[,ia,drop=F]
  n <- length(ia)
  sa <- integer(n) # sequence of states
  na <- 0 # counter of states
  for (j in 1:n) {
    if (sa[j]!=0) next
    sa[j] <- (na <- na+1)
    if (j<n) for (k in (j+1):n) {
      if (all(v[,j]==v[,k])) sa[k] <- na
    }
  }
  a <- matrix(F,nrow(sm$s),na) # the states (available species)
  rownames(a) <- rownames(sm$s)
  p <- a # the states (potentially available species... to be corrected using closure)
  for (j in 1:na) {
    k <- match(j,sa) # first occurrence of state j in sa
    a[,j] <- v[,k]
    p[rn.closure(rn,which(v[,k])),j] <- T # the closure of the state
  }
  ap <- integer(ncol(p)) # abstract states corresponding closed state
  np <- 0 # counter of states
  for (j in 1:ncol(p)) {
    if (ap[j]!=0) next
    ap[j] <- (np <- np+1)
    if (j<ncol(p)) for (k in (j+1):ncol(p)) {
      if (all(p[,j]==p[,k])) ap[k] <- np
    }
  }
  j <- match(1:np,ap)
  p <- p[,j,drop=F]
  sp <- ap[sa]
  if (length(sp)>1) j <- c(1,1+which((sp[-1]!=sp[-length(sp)])))
  ip <- ia[j]
  sp <- sp[j]
  return(list(ia=ia,sa=sa,a=a,ap=ap,ip=ip,sp=sp,p=p))
}

# returns tconv (convergence time) an ttot (total time) for simulation results sm
sm.conv <- function(sm) {
  f <- sign(sm$v[,ncol(sm$v)])
  m <- apply(sign(sm$v),2,function(v) all(v==f))
  cnv <<- m
  k <- match(F,rev(m))
  return(c(tconv = if (!is.na(k)) sm$t[length(sm$t)-k+1] else 0, ttot = sm$t[length(sm$t)]))
}

# returns the final closed state of reaction network rn and simulation results sm
sm.final <- function(rn,sm) {
  v <- sm$v[,ncol(sm$v)]
  return(rn.support(rn,which(v>0)))
}

sm.display <- function(sm,L=1000) {
  require(ggplot2)
  require(reshape2)
  require(cowplot)
  l <- length(sm$t)
  if (l>2*L) {
    s <- floor(l/L)
    i <- seq(1,l,s)
    if (i[length(i)]!=l) i <- c(i,l)
  }
  else
    i <- 1:l
  d.s <- melt(t(sm$s[,i]),varnames=c("time","species"))
  d.s$time <- sm$t[i[d.s$time]]
  d.v <- melt(t(sm$v[,i]),varnames=c("time","reaction"))
  d.v$time <- sm$t[i[d.v$time]]
  g1 <- ggplot(d.s, aes(x=time, y=value, col=species)) + geom_line()
  g2 <- ggplot(d.v, aes(x=time, y=value, col=reaction)) + geom_line()
  plot_grid(g1,g2,ncol=1,nrow=2)
}

sm.example <- function(rn=sm.genrn(12),n=1000,dt=.1) {
  cat(nrow(rn$mr),"species,",ncol(rn$mr),"reactions:\n")
  rn.display(rn)
  o <- rn.linp_org(rn)
  cat("needed inflow",o$ifl,"\n")
  cat("overproducible",o$ovp,"\n")
  sm <<- sm.sim(rn,n=n,dt=dt,t0=0,momentum=0,norm=0)
  print(sm.conv(sm))
  print(sm$s[,ncol(sm$s)])
  print(sm$v[,ncol(sm$v)])
  print(floor(log10(pmax(0,sm$s[,ncol(sm$s)]))))
  print(floor(log10(sm$v[,ncol(sm$v)])))
  if (!is.null(rn$sid)) cat("final state:",rn$sid[sm.final(rn,sm)],"\n")
  else cat("final state:",sm.final(rn,sm),"\n")
  gc()
  sm.display(sm)
}
