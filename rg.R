
# random generation of reaction networks

# shift left a vector preserving the last element
rg.shl <- function(v) c(v[-1],v[length(v)])

# returns a random reaction template
# up to N[1] reactants and N[2] products (depending on pnot),
# pcat: catalyzer prob., pmix: mixed reac/prod p., pmp: produced mixed p., ppr: pure reac. p., ppp: pure prod. p.
# There is at least one consumed and one produced species.
# returned structure: list(nc,nm,np) (# catalyzers, # mixed consumed/produced species, # pure consumed/produced species)
rg.rrt <- function(N=c(2,2),pcat=c(0.1,0),pmix=0.1,pmp=0.5,ppr=0.8,ppp=0.8) {
  r <- list(nc=0,nm=c(0,0),np=c(0,0))
  # catalyzers species
  i <- min(N)
  while(pcat[1]>0 && i>0 && all(N>1)) {
    if (runif(1)<pcat[1]) { r$nc <- r$nc+1; N <- N - 1 }
    pcat <- rg.shl(pcat)
    i <- i - 1
  }
  # mixed species
  i <- min(N) - (sum(N)==2) # because at least one consumed and one produced species
  while(pmix[1]>0 && i>0 && all(N>0)) {
    if (runif(1)<pmix[1]) {
      if (runif(1)<pmp) { # produced
        if (N[1]>1 || r$nm[1]>0) { r$nm <- r$nm + c(0,1); N <- N - 1 }
      } else { # consumed
        if (N[2]>1 || r$nm[2]>0) { r$nm <- r$nm + c(1,0); N <- N - 1 }
      }
    }
    pmix <- rg.shl(pmix)
    i <- i - 1
  }
  r$np <- as.numeric(r$nm==0)  # to ensure at least one consumed and one produced species
  if (r$np[1]==1) ppr <- rg.shl(ppr)  # one reactant added
  if (r$np[2]==1) ppp <- rg.shl(ppp)  # one product added
  N <- N - r$np
  # remaining pure species
  if (ppp==0) N[2] <- 0  # no more product species allowed
  else if (ppp==1) N[1] <- 0  # no more reactant species allowed
  # pure reactants
  i <- N[1]
  while(ppr[1]>0 && i>0) {
    if (runif(1)<ppr[1]) { N[1] <- N[1] - 1; r$np[1] <- r$np[1] + 1 }
    ppr <- rg.shl(ppr)
    i <- i - 1
  }
  # pure products
  i <- N[2]
  while(ppp[1]>0 && i>0) {
    if (runif(1)<ppp[1]) { N[2] <- N[2] - 1; r$np[2] <- r$np[2] + 1 }
    ppp <- rg.shl(ppp)
    i <- i - 1
  }
  return(r)
}

# returns a reaction template selected randomly from a given weighed set (each element is the row of a matrix)
# columns: catalyzers, mixed consumed, mixed produced, pure consumed, pure produced, weight 
rg.rsrt <- function(m=rbind( c(1,0,0,1,1,10),
                             c(1,0,1,1,0,10),
                             c(0,0,0,1,1,80)
                            )) {
  p <- m[,6]; p <- cumsum(p/sum(p)) # cumulative probability
  i <- match(T,runif(1)<p)
  return(list(nc=m[i,1],nm=m[i,2:3],np=m[i,4:5]))
}

# returns a data frame with all cardinality-mass combinations and joint distribution:
#   species cardinality C, masses M, cardinality distribution dC, mass distribution dM
rg.dCM <- function(C=1:4, M=1:4, dC=rep(1,length(C)), dM=rep(1,length(M))) {
  df <- merge(cbind(C=C),cbind(M=M))
  df$d <- abs(dC[df$C]*dM[df$M])
  df$d <- df$d/sum(df$d)
  return(df)
}

# generates a generating function for random reactions
# dCM: dataframe of cardinality-mass combinations, frt: function for creating reaction templates
# The result is a function whose result is a matrix were each column is left card., right card. and mass for a species
rg.genreac <- function(dCM=rg.dCM(), frt=rg.rrt) {
  
  # selection of a cardinality/mass combination, any cardinality but n, only m mass, only mt total mass (card*mass)  
  select <- function(n=NULL,m=NULL,mt=NULL) {
    d <- dCM$d
    if (!is.null(n)) d[dCM$C==n] <- 0
    if (!is.null(m)) d[dCM$M!=m] <- 0
    if (!is.null(mt)) d[dCM$C*dCM$M!=mt] <- 0
    s <- sum(d>0)
    if (s==0) return(NULL)
    if (s==1) return(dCM[which(d>0),1:2])
    return(dCM[match(T,runif(1)<cumsum(d/sum(d))),1:2])
  }
  
  rg.unbalreac <- function(rt) {
    r <- matrix(0,3,sum(rt$nc,rt$nm,rt$np))
    if (rt$nc>0) for (i in 1:rt$nc) { # catalyzers
      s <- select()
      r[,i] <- c(s$C,s$C,s$M)
    }
    j <- rt$nc
    if (rt$nm[1]>0) for (i in 1:rt$nm[1]) { # mixed consumed species
      repeat { s1 <- select(); s2 <- select(n=s1$C,m=s1$M); if (!is.null(s2)) break }
      if (s1$C>s2$C) r[,i+j] <- c(s1$C,s2$C,s1$M) else r[,i+j] <- c(s2$C,s1$C,s1$M)
    }
    j <- j + rt$nm[1]
    if (rt$nm[2]>0) for (i in 1:rt$nm[2]) { # mixed produced species
      repeat { s1 <- select(); s2 <- select(n=s1$C,m=s1$M); if (!is.null(s2)) break }
      if (s1$C<s2$C) r[,i+j] <- c(s1$C,s2$C,s1$M) else r[,i+j] <- c(s2$C,s1$C,s1$M)
    }
    j <- j + rt$nm[2]
    if (rt$np[1]>0) for (i in 1:rt$np[1]) { # pure reactants
      s <- select()
      r[,i+j] <- c(s$C,0,s$M)
    }
    j <- j + rt$np[1]
    if (rt$np[2]>0) for (i in 1:rt$np[2]) { # pure products
      s <- select()
      r[,i+j] <- c(0,s$C,s$M)
    }
    return(r)
  }

  rg.reac <<- function(rt=frt(),balanced=T) {
    r <- rg.unbalreac(rt)
    if (!balanced) return(r)
    b <- (r[2,]-r[1,])*r[3,]
    B <- sum(b)
    if (B!=0) repeat {
      ra <- rg.unbalreac(rt)
      ba <- (ra[2,]-ra[1,])*ra[3,]
      Ba <- sum(ba)
      i <- which.min(abs(B+c(0,ba-b)))
      ia <- which.min(abs(Ba+c(0,b-ba)))
      if (i>1) B <- B + ba[i-1] - b[i-1]
      if (ia>1) Ba <- Ba + b[ia-1] - ba[ia-1]
      if (abs(B)<abs(Ba)) {
        if (i>1) { b[i-1] <- ba[i-1]; r[,i-1] <- ra[,i-1] }
      }
      else {
        if (ia>1) { ba[ia-1] <- b[ia-1]; ra[,ia-1] <- r[,ia-1] }
        r <- ra; b <- ba; B <- Ba
      }
      if (sum(B)==0) break
    }
    return(r)
  }

  return(invisible(rg.reac))
}

# returns a list of L reactions with U unbalanced ones using the generating function g
rg.lreac <- function(L=100,U=0,g=rg.genreac()) {
  lapply(1:L,function(i) g(balanced=(i>U)))
}

# returns an under construction reaction networks object to be initialized by a list of reactions 
rg.ucnet <- setRefClass("ucnet",fields=c("mr","mp","w","u","scl","s","rs"))
rg.ucnet$methods(
  initialize = function(lreac=NULL) { # called when creating the object 'rg.ucnet(...)'
    if (is.null(lreac)) return(0)
    Ns <- sum(sapply(lreac,function(r) ncol(r))) # number of species
    mr <<- mp <<- matrix(0,Ns,length(lreac)) # reactant an product matrices
    w <<- numeric(Ns) # species weight
    u <<- matrix(0,Ns,3) # how many times the species is used as (1: consumed, 2: produced, 3: catalyst)
    s <<- 1:Ns # species substitution, initially every species is just itself
    rs <<- 1:Ns # remaining species, initially all of them
    l <- numeric(Ns) # species reaction index
    i <- 0; j <- 1
    for (r in lreac) {
      k <- i + (1:ncol(r))
      mr[k,j] <<- r[1,]; mp[k,j] <<- r[2,]
      u[k,1] <<- r[1,]>r[2,]; u[k,2] <<- r[1,]<r[2,]; u[k,3] <<- r[1,]==r[2,]
      w[k] <<- r[3,]; l[k] <- j
      i <- i + ncol(r); j <- j + 1
    }
    scl <<- lapply(1:Ns, function(i) which(l[i]!=l & w[i]==w)) # compatibles species for each species
  },

  unify = function(i,j) { # tries to unify species i and j, returns T if successful F if not
    if (i<1 || j<1 || i>length(w) || j> length(w) || !(j %in% scl[[i]])) return(F)
    if (i>j) { k <- i; i <- j; j <- k } # always i < j
    I.i <- which(mr[i,]+mp[i,]>0) # involved reactions (i species)
    I.j <- which(mr[j,]+mp[j,]>0) # involved reactions (j species)
    ok <- T
    for (ii in I.i) {
      k <- which(mr[j,I.j]==mr[i,ii] & mp[j,I.j]==mp[i,ii])
      for (jj in I.j[k])
        if (all(mr[-c(i,j),ii]==mr[-c(i,j),jj]&&mp[-c(i,j),ii]==mp[-c(i,j),jj])) { ok <- F; break }
      if (!ok) break
    }
    if (ok) { # the unification is ok: it doesn't create a redundant reaction
      s[rs[j]] <<- rs[i] # i and j are unified by substituting j by i
      rs <<- rs[-j]
      mr[i,] <<- mr[i,] + mr[j,]; mr <<- mr[-j,]
      mp[i,] <<- mp[i,] + mp[j,]; mp <<- mp[-j,]
      w <<- w[-j]
      u[i,] <<- u[i,] + u[j,]; u <<- u[-j,]
      scl[[i]] <<- intersect(scl[[i]],scl[[j]]); scl[[j]] <<- NULL
      scl <<- lapply( scl, function(sc) {
                             if (i %in% sc) if (!j %in% sc) sc <- setdiff(sc,i)
                             sc <- setdiff(sc,j)
                             sc <- sc-(sc>j)
                             return(sc)
                           } )
    }
    else { scl[[i]] <<- setdiff(scl[[i]],j); scl[[j]] <<- setdiff(scl[[j]],i) } # not ok: the compatibility is erased
    return(ok)
  },
  
  rn = function() list(w=w,mr=mr,mp=mp,s=s)
)

# attempts to cyclically connect a reaction network under construction to form a semi-organization
# connections are achieved through the unification of species
rg.cyclcon <- function(ucn) {
  u <- which(ucn$u[,1]>0 & ucn$u[,2]==0 | ucn$u[,1]==0 & ucn$u[,2]>0) # consumed-only or produced-only species
  while(length(u)>0) {
    i <- -sample(-u,1)
    cs <- ucn$scl[[i]]
    cs <- cs[ucn$u[cs,1+(ucn$u[i,1]>0)]>0]
    if (length(cs)==0) { u <- setdiff(u,i); next } # no unification available for i
    j <- cs[ucn$u[cs,1+(ucn$u[i,2]>0)]==0]
    if (length(j)>0) j <- -sample(-j,1) else j <- -sample(-cs,1)
    if (ucn$unify(i,j)) {
      u <- setdiff(u,c(i,j))
      if (length(u)==0) break
      u <- u - (u>max(i,j))
    }
  }
}

# makes N new random unifications to connect a reaction network under construction
rg.rcon <- function(ucn,N=1) {
  while (N>0) {
    s <- which(sapply(ucn$scl,function(e) length(e)>0))
    if (length(s)==0) break
    i <- -sample(-s,1)
    j <- -sample(-ucn$scl[[i]],1)
    if (ucn$unify(i,j)) N <- N-1
  }
}

rg.example <- function(L=10,s=1) {
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
  rn.display(ucn$rn())
  o <- rn.linp_org(ucn$rn())
  cat("needed inflow",o$ifl,"\n")
  cat("overproducible",o$ovp,"\n")
  "The result is in 'ucn' global variable"
}
