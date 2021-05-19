
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
# minimum activation threshold e (vector: low, high) and maximum species concentration w
sm.sim <- function(rn,n=1000,dt=.1,s0=runif(nrow(rn$mr),.9,1.1),t0=0,p=runif(ncol(rn$mr),.9,1.1),e=1e-3*c(1,10), w=1e2) {
  s <- matrix(0,nrow(rn$mr),n)  # state matrix (species concentration, one column vector per iteration)
  if (!is.null(rn$sid)) row.names(s) <- rn$sid else row.names(s) <- paste0("s",1:nrow(rn$mr))
  v <- matrix(0,ncol(rn$mr),n) # process matrix (reactions activity, one column vector per iteration)
  row.names(v) <- paste0("r",1:ncol(rn$mr))
  t <- rep(0,n) # time
  m <- rn$mp - rn$mr
  if (length(e)<2) e <- c(e,e)
  if (t0<=0) s[,1] <- pmax(0,s0)
  cs <- logical(nrow(rn$mr)) # current state (species available)
  for (i in 1:n) {
    if (i>1) {
      ds <- m %*% v[,i-1,drop=F] # change of concentrations per time unit
      k <- which(ds<0) # species that could drop to 0
      f <- min(-s[k,i-1]/ds[k],1)
      t[i] <- t[i-1] + f*dt
      s[,i] <- pmin(pmax(0,s[,i-1]+f*ds*dt),w)
    }
    cs[s[,i]>e[2]] <- T # species over high threshold are available
    cs[s[,i]<=e[1]] <- F # species below low threshold are not available
    k.s <- which(cs) # available species
    k.r <- which(rbind(cs==F) %*% rn$mr == 0) # reactions with all reactants available
    if (length(k.r)==0) next
    v[k.r,i] <- dt*p[k.r]*exp(rbind(log(s[k.s,i])) %*% rn$mr[k.s,k.r])
  }
  return(list(s=s,v=v,t=t))
}

# returns tconv (convergence time) an ttot (total time) for simulation results sm
sm.conv <- function(sm) {
  f <- sign(sm$v[,ncol(sm$v)])
  m <- apply(sign(sm$v),2,function(v) all(v==f))
  cnv <<- m
  k <- match(F,rev(m))
  return(c(tconv = if (!is.na(k)) sm$t[length(sm$t)-k+1] else 0, ttot = sm$t[length(sm$t)]))
}

# returns the final state of reaction network rn and simulation results sm
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

sm.example <- function(rn=sm.genrn(),n=1000,dt=.1) {
  cat(nrow(rn$mr),"species,",ncol(rn$mr),"reactions:\n")
  rn.display(rn)
  o <- rn.linp_org(ucn$rn())
  cat("needed inflow",o$ifl,"\n")
  cat("overproducible",o$ovp,"\n")
  sm <<- sm.sim(rn,n=n,dt=dt)
  print(sm.conv(sm))
  print(sm$s[,ncol(sm$s)])
  print(sm$v[,ncol(sm$v)])
  print(floor(log10(sm$s[,ncol(sm$s)])))
  print(floor(log10(sm$v[,ncol(sm$v)])))
  if (!is.null(rn$sid)) cat("final state:",rn$sid[sm.final(rn,sm)],"\n")
  else cat("final state:",sm.final(rn,sm),"\n")
  gc()
  sm.display(sm)
}
