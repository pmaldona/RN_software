
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

# Simulation of a reaction network rn for n iterations with maximal random reactions.
# Initial state s0 (vector, one component per species), reactions preference p (one component per reaction),
sm.mrrsim <- function( rn, n=100*ncol(rn$mr), s0=runif(nrow(rn$mr),.9,1.1), p=runif(ncol(rn$mr),.3,3) ) {
  s <- matrix(0,nrow(rn$mr),n) # state matrix (species concentration, one column vector per iteration)
  sa <- matrix(F,nrow(rn$mr),n) # state matrix abstraction
  rownames(s) <- rownames(rn$mr)
  v <- matrix(0,ncol(rn$mr),n) # process matrix (reactions activity, one column vector per iteration)
  rownames(v) <- colnames(rn$mr)
  t <- numeric(n) # time vector
  m <- rn$mp - rn$mr
  cs <- logical(nrow(rn$mr)) # current state (boolean vector of the species that are available)
  s[,1] <- s0 # starting state
  for (i in 1:n) {
    if (i>1) {
      k <- sample(nrow(v),1,prob=v[,i-1])  # randomly chosen reaction
      l <- which(m[,k]<0)  # species that are consumed
      if (length(l)>0) f <- min(-s[l,i-1]/m[l,k])  # how much the reaction will be applied
      else f <- 1  # just once if no consumption
      t[i] <- t[i-1] + 1/nrow(v)
      s[,i] <- s[,i-1]+f*m[,k]
    }
    cs[] <- s[,i] > 0
    sa[,i] <- cs
    k.r <- which(rbind(cs==F) %*% rn$mr == 0) # reactions with all reactants available
    if (length(k.r)==0) break
    v[k.r,i] <- p[k.r]
  }
  return(list(s=s,sa=sa,v=v,t=t))  # s: species concentration, sa: species abstract existence, v: flow vector, t: time
}

# Simple simulation of a reaction network rn for n iterations at dt time step with a constant flow vector.
# Initial state s0 (vector, one component per species), flow vector p (one component per reaction),
# p may be normalized, deficit of species allowed.
sm.cfsim <- function( rn, n=1000, dt=.1, s0=runif(nrow(rn$mr),.9,1.1), p=runif(ncol(rn$mr),.3,3), norm=F ) {
  if (norm) p <- p/mean(p)
  s <- matrix(0,nrow(rn$mr),n) # state matrix (species concentration, one column vector per iteration)
  sa <- matrix(F,nrow(rn$mr),n) # state matrix abstraction
  rownames(s) <- rownames(rn$mr)
  v <- matrix(0,ncol(rn$mr),n) # process matrix (reactions activity, one column vector per iteration)
  rownames(v) <- colnames(rn$mr)
  t <- numeric(n) # time vector
  m <- rn$mp - rn$mr
  cs <- logical(nrow(rn$mr)) # current state (boolean vector of the species that are available)
  s[,1] <- s0 # starting state
  for (i in 1:n) {
    if (i>1) {
      ds <- m %*% v[,i-1] # change of concentrations per time unit
      t[i] <- t[i-1] + dt
      s[,i] <- s[,i-1]+ds*dt
    }
    cs[] <- s[,i] > 0
    sa[,i] <- cs
    k.s <- which(cs) # available species
    k.r <- which(rbind(cs==F) %*% rn$mr == 0) # reactions with all reactants available
    if (length(k.r)==0) next
    v[k.r,i] <- p[k.r]
  }
  return(list(s=s,sa=sa,v=v,t=t))  # s: species concentration, sa: species abstract existence, v: flow vector, t: time
}

# Euler method simulation of a mass action kinetics of a reaction network rn for n iterations at dt time step
# with initial state s0 (vector of concentrations, one component per species), mass action parameter p (vector,
# one component per reaction), minimum concentration threshold e, .
# Inflow species are kept constant if inflow parameter is true.
# The time step is strictly respected so that concentrations can eventually drop below zero.
# The time unit is warped so that the maximum concentration drop is 1 per time unit.
# The dt time step can be interpreted as the maximum concentration drop per iteration.
sm.maksim <- function(rn,n=1000,dt=.2,s0=runif(nrow(rn$mr),.9,1.1),p=runif(ncol(rn$mr),.9,1.1),e=.1,a=1,inflow=T) {
  s <- matrix(0,nrow(rn$mr),n) # state matrix (species concentration, one column vector per iteration)
  s[,1] <- s0 # initial state
  sa <- matrix(F,nrow(rn$mr),n) # state matrix abstraction
  rownames(s) <- rownames(rn$mr)
  rownames(sa) <- rownames(rn$mr)
  v <- matrix(0,ncol(rn$mr),n) # process matrix (reactions activity, one column vector per iteration)
  rownames(v) <- colnames(rn$mr)
  t <- numeric(n) # time vector
  m <- rn$mp - rn$mr # process matrix
  cs <- logical(nrow(rn$mr)) # current state (species available)
  if (inflow==F)
    inflow <- numeric(0)
  else if (inflow==T)
    inflow <- which(rowSums(rn$mp[,which(colSums(rn$mr)==0),drop=F])>0) # else a given vector of species (index)
  for (i in 1:n) {
    if (i>1) {
      ds <- m %*% v[,i-1] # change of concentrations per normalized time unit
      mds <- max(-ds,ds/10) # normalization factor for ds (minimum component will be set to -1, maximum to 10)
      if (mds==0) mds <- 1
      t[i] <- t[i-1] + dt*exp(-mlv)/mds
      s[,i] <- s[,i-1]+ds*dt/mds
    }
    cs[] <- s[,i] > e
    sa[,i] <- cs
    k.s <- which(cs) # currently available species
    k.r <- which(rbind(cs==F) %*% rn$mr == 0 & p>0) # reactions with all reactants available and positive constant
    if (length(k.r)==0) break  # we have reached a non reactive end state, nothing more to simulate...
    if (a==0) lv <- log(p[k.r])
    else lv <- log(p[k.r]) + a*rbind(log(s[k.s,i])) %*% rn$mr[k.s,k.r]  # the flow vector in logarithmic scale
    mlv <- max(lv)  # normalization factor for flow vector and time unit
    v[k.r,i] <- exp(lv-mlv)  # normalized flow vector (maximum component is 1)
  }
  return(list(s=s[,1:i,drop=F],sa=sa[,1:i,drop=F],v=v[,1:i,drop=F],t=t[1:i]))
  # s: species concentration, sa: species abstract existence, v: flow vector, t: time
}

# Euler method simulation of a reaction network rn for n iterations at dt time step (adaptative)
# with initial state s0 (infused up to time t0), mass action parameter p (vector, one component per reaction),
# minimum activation threshold e (vector: low, high), maximum species concentration w, acceptance of deficit of species,
# momentum (0 => no momentum), normalization norm (norm==0 => no normalization), inflow species kept constant.
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
  k <- match(F,rev(m))
  return(c(tconv = if (!is.na(k)) sm$t[length(sm$t)-k+1] else 0, ttot = sm$t[length(sm$t)]))
}

# returns the final closed state (species) of reaction network rn and simulation results sm
sm.final <- function(rn,sm) {
  v <- sm$v[,ncol(sm$v)]
  s <- rn.support(rn,which(v>0))
  return(rn.closure(rn,s))
}

sm.display <- function(sm,L=1000,simple=F) {
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
  if (simple) {
    d <- cbind(species=colSums(sm$sa[,i]),reactions=colSums(sm$v[,i]>0))
    d <- melt(d,varnames=c("time","number of"))
    d$iteration <- i
    d <<- d
    ggplot(d,aes(x=iteration, y=value, col=`number of`)) + geom_line()
  }
  else {
    d.s <- melt(t(sm$s[,i]),varnames=c("time","species"))
    d.v <- melt(t(sm$v[,i]),varnames=c("time","reaction"))
    d.s$time <- sm$t[i[d.s$time]]
    d.v$time <- sm$t[i[d.v$time]]
    g1 <- ggplot(d.s, aes(x=time, y=value, col=species)) + geom_line()
    g2 <- ggplot(d.v, aes(x=time, y=value, col=reaction)) + geom_line()
    plot_grid(g1,g2,ncol=1,nrow=2)
  }
}

sm.example <- function(rn=sm.genrn(12),n=1000,dt=.1,e=.2,a=1) {
  cat(nrow(rn$mr),"species,",ncol(rn$mr),"reactions:\n")
  rn.display(rn)
  o <- rn.linp_org(rn)
  cat("needed inflow",o$ifl,"\n")
  cat("overproducible",o$ovp,"\n")
  sm <<- sm.maksim(rn,n=n,dt=dt,e=e,a=a)
  # sm <<- sm.sim(rn,n=n,dt=dt,t0=0,e=c(.2,.2),w=Inf,momentum=0)
  # sm <<- sm.sim(rn,n=n,dt=dt,t0=0,momentum=0)
  # sm <<- sm.cfsim(rn,n=n,dt=dt)
  # sm <<- sm.mrrsim(rn,n=n)
  print(sm.conv(sm))
  print(sm$s[,ncol(sm$s)])
  print(sm$v[,ncol(sm$v)])
  print(floor(log10(pmax(0,sm$s[,ncol(sm$s)]))))
  print(floor(log10(sm$v[,ncol(sm$v)])))
  if (!is.null(rn$sid)) cat("final state:",rn$sid[sm.final(rn,sm)],"\n")
  fs <- sm.final(rn,sm)
  cat("final species state:", if (length(fs)>0) fs else "void","\n")
  if (length(fs)>0) {
    rnf <- rn.sub(rn,fs)
    cat(nrow(rnf$mr),"species,",ncol(rnf$mr),"reactions:\n")
    rn.display(rnf)
    rnf <<- rnf
  }
  gc()
  g <- sm.display(sm,simple=T)
  print(g)
  Sys.sleep(1)
}
