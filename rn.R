
# reaction networks
#
# representation: list(sid,w,mr,mp)  (species identifier (name), species weight, matrix of reactants, matrix of products)
# additional components (sx,rx)  (species index, identifier index) needed for mapping sub-networks

library("limSolve")

# generate simple reaction networks to make tests
rn.testrn <- function(i=1) {
  if (i==1) { # farm model
    sid <- c("water","grass","cows","infrastructure","milk","dung","worms",
             "fertilizer","chickens","eggs","grain","straw","farmer","money")
    mr <- matrix(0,length(sid),17); rownames(mr) <- sid; colnames(mr) <- paste0("r",1:17); mp <- mr
    mp["water",1] <- 1
    mr[c("grass","cows","infrastructure","water"),2] <- 1; mp[c("milk","cows","dung","infrastructure"),2] <- 1
    mr["dung",3] <- 1; mp[c("worms","fertilizer"),3] <- 1
    mr[c("grass","worms","chickens","infrastructure"),4] <- 1; mp[c("chickens","eggs","fertilizer","infrastructure"),4] <- 1
    mr[c("grain","chickens","infrastructure"),5] <- 1; mp[c("chickens","eggs","fertilizer","infrastructure"),5] <- 1
    mr[c("water","fertilizer","grass"),6] <- 1; mp["grass",6] <- 2
    mr[c("eggs","grass","worms","infrastructure"),7] <- 1; mp[c("chickens","fertilizer","infrastructure"),7] <- 1
    mr[c("fertilizer","water","grain","infrastructure"),8] <- 1; mp[c("grain","straw","infrastructure"),8] <- c(10,1,1)
    mr[c("straw","cows","water","infrastructure"),9] <- 1; mp[c("milk","cows","infrastructure"),9] <- 1
    mr[c("eggs","farmer"),10] <- 1; mp[c("money","farmer"),10] <- 1
    mr[c("milk","farmer"),11] <- 1; mp[c("money","farmer"),11] <- 1
    mr[c("grain","farmer"),12] <- 1; mp[c("money","farmer"),12] <- 1
    mr[c("straw","farmer"),13] <- 1; mp[c("money","farmer"),13] <- 1
    mr[c("chickens","farmer"),14] <- 1; mp[c("money","farmer"),14] <- 1
    mr["worms",15] <- 1
    mr["infrastructure",16] <- 1
    mr[c("money","farmer","infrastructure"),17] <- 1; mp[c("farmer","infrastructure"),17] <- c(1,2)
    return(list(mr=mr,mp=mp))
  }
  else if (i==2) {
    sid <- c("a","b","c")
    mr <- cbind(c(1,0,0),c(0,1,0),c(0,0,1),c(1,0,0))
    mp <- cbind(c(0,1,0),c(0,0,1),c(1,0,0),c(0,0,1))
    rownames(mr) <- rownames(mp) <- c("a","b","c")
    colnames(mr) <- colnames(mp) <- paste0("r",1:4)
    return(list(mr=mr,mp=mp))
  }
  else if (i==3) { # SEIRS model for contagion
    sid <- c("s","e","i","r")
    mr <- cbind(c(1,0,1,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1))
    mp <- cbind(c(0,1,1,0),c(0,0,1,0),c(0,0,0,1),c(1,0,0,0))
    rownames(mr) <- rownames(mp) <- c("s","e","i","r")
    colnames(mr) <- colnames(mp) <- paste0("r",1:4)
    return(list(mr=mr,mp=mp))
  }
  else if (i==4) { # mini ecology
    sid <- c("s","p","h","c") # soil, plant, herbivore, carnivore
    mr <- cbind(c(1,1,0,0),c(0,1,0,0),c(0,1,1,0),c(0,0,1,0),c(0,0,1,1),c(0,0,0,1))
    mp <- cbind(c(0,2,0,0),c(1,0,0,0),c(1,0,2,0),c(0,0,0,0),c(0,0,0,2),c(0,0,0,0))
    rownames(mr) <- rownames(mp) <- c("s","p","h","c") # soil, plant, herbivore, carnivore
    colnames(mr) <- colnames(mp) <- paste0("r",1:6)
    return(list(mr=mr,mp=mp))
  }
}

# read a reaction network expressed in Antimony text format
rn.read <- function(file) {
  rematch <- function(s,p) { 
    i <- regexpr(p,s)
    if (i==-1) return("")
    else return(substr(s,i,i+attr(i,"match.length")-1))
  }
  l <- readLines(file)
  lr <- list(); lp <- list()
  for (s in l) {
    s <- gsub("[ \t]","",s)
    s <- gsub("^.*:","",s)
    s <- gsub("[;#].*$","",s)
    op <- rematch("->|=>",s)
    s <- strsplit(s,"->|=>")[[1]]
    if (length(s)!=2) next
    s1 <- strsplit(s[1],"+",fixed=T)[[1]]
    s2 <- strsplit(s[2],"+",fixed=T)[[1]]
    if (length(s1)==0)
      er <- list()
    else
      er <- list( c = sapply(sapply(s1,rematch,"[0-9]*"), function(s) if(s=="") 1 else as.integer(s)),
                  id = sapply(s1,rematch,"[^0-9].*") )
    if (length(s2)==0)
      ep <- list()
    else
      ep <- list( c = sapply(sapply(s2,rematch,"[0-9]*"), function(s) if(s=="") 1 else as.integer(s)),
                  id = sapply(s2,rematch,"[^0-9].*") )
    lr <- c(lr,list(er))
    lp <- c(lp,list(ep))
    if (op=="->") {
      lr <- c(lr,list(ep))
      lp <- c(lp,list(er))
    }
  }
  sid <- c()
  for (i in 1:length(lr)) sid <- union(sid,union(lr[[i]]$id,lp[[i]]$id))
  mr <- mp <- matrix(0,length(sid),length(lr))
  rownames(mr) <- rownames(mp) <- sid
  colnames(mr) <- colnames(mp) <- paste0("r",1:ncol(mr))
  for (i in 1:length(lr)) {
    if (length(lr[[i]])>0) mr[sapply(lr[[i]]$id,match,sid),i] <- lr[[i]]$c
    if (length(lp[[i]])>0) mp[sapply(lp[[i]]$id,match,sid),i] <- lp[[i]]$c
  }
  list(mr=mr,mp=mp)
}

# merges the reactions of one or more reaction networks 
rn.merge <- function(...) {
  mr <- NULL; mp <- NULL
  if (...length()>0) for (i in 1:...length()) {
    mr1 <- mr; mp1 <- mp; sid1 <- rownames(mr)
    nc1 <- ncol(mr1); if (is.null(nc1)) nc1 <- 0
    mr2 <- ...elt(i)$mr; mp2 <- ...elt(i)$mp; sid2 <- rownames(mr2)
    nc2 <- ncol(mr2); if (is.null(nc2)) nc2 <- 0
    sid <- sort(union(sid1,sid2))
    mr <- matrix(0,length(sid),nc1+nc2)
    rownames(mr) <- sid
    mp <- mr
    if (nc1>0) {
      mr[sid1,1:nc1] <- mr1
      mp[sid1,1:nc1] <- mp1
    }
    if (nc2>0) {
      mr[sid2,1:nc2+nc1] <- mr2
      mp[sid2,1:nc2+nc1] <- mp2
    }
  }
  i <- which(rowSums(mr)+rowSums(mp)==0)
  if (length(i)>0) { mr <- mr[-i,,drop=F]; mp <- mp[-i,,drop=F] }  # unused species are eliminated
  i <- which(sapply(1:ncol(mr), function(k) all(mr[,k]==mp[,k])))
  if (length(i)>0) { mr <- mr[,-i,drop=F]; mp <- mp[,-i,drop=F] }  # null reactions are eliminated
  i <- which(duplicated(mr,MARGIN=2) & duplicated(mp,MARGIN=2))
  if (length(i)>0) { mr <- mr[,-i,drop=F]; mp <- mp[,-i,drop=F] }  # duplicated reactions are eliminated
  if (ncol(mr)>0) {
    colnames(mr) <- colnames(mp) <- paste0("r",1:ncol(mr))
  }
  return(list(mr=mr,mp=mp))
}

# displays selected reactions i (defaults to all) of a reaction network rn
# if file is not "" instead of displaying the network it is saved in that file
rn.display <- function(rn,i=1:ncol(rn$mr),file="") {
  if (file!="" && file.exists(file)) file.remove(file)
  if (!is.null(rownames(rn$mr))) sid <- rownames(rn$mr)
  else sid <- paste0("s",1:nrow(rn$mr))
  if (!is.null(colnames(rn$mr))) rid <- colnames(rn$mr)
  else rid <- paste0("r",1:ncol(rn$mr))
  if (length(i)>0) for (k in i) {
    k.r <- which(rn$mr[,k]>0)
    k.p <- which(rn$mp[,k]>0)
    m <- cbind(rbind(rn$mr[k.r,k],sid[k.r]),rbind(rn$mp[k.p,k],sid[k.p]))
    m <- rbind(m," + ")
    if(length(k.r)>0) m[3,length(k.r)] <- " => "
    else m <- cbind(c("","∅"," => "),m)
    if(length(k.p)==0) m <- cbind(m,c("","∅",""))
    m[3,ncol(m)] <- ";\n"
    m[1,m[1,]=="1"] <- ""
    cat(rid[k],":\t",m,sep="",file=file,append=T)
  }
}

# returns the species of a reaction network rn supporting a set of reactions r (defaults to all reactions)
rn.support <- function(rn,r=1:ncol(rn$mr)) return(which(rowSums(rn$mr[,r,drop=F]+rn$mp[,r,drop=F])>0))

# returns the reactions of rn that are supported by a set of species s (defaults to all species)
rn.supported <- function(rn,s=1:nrow(rn$mr)) return(which(colSums(rn$mr[s,,drop=F]>0)==colSums(rn$mr>0)))

# returns a sub-network of a reaction network nr by species s (defaults to all), unreactive species are trimmed
rn.sub <- function(rn,s=1:nrow(rn$mr)) {
  m <- rn$mr + rn$mp # matrix of species used by reactions
  rs <- which(colSums(m)==colSums(m[s,,drop=F])) # reactions supported by the selected species
  sr <- which(rowSums(m[,rs,drop=F])>0) # species used in the selected reactions
  srn <- list()
  if (!is.null(rn$sid)) srn$sid <- rn$sid[sr]
  if (!is.null(rn$w)) srn$w <- rn$w[sr]
  srn$mr <- rn$mr[sr,rs]
  srn$mp <- rn$mp[sr,rs]
  return(srn)
}

# returns independently transient species
rn.trans <- function(rn,mr=rn$mr,mp=rn$mp) {
  i <- which(colSums(mr>0)==1)  # reactions with only one reactant
  return(which(rowSums(mr[,i,drop=F])>0 & rowSums(mp-mr>0)==0))  # lone reactants that are not produced
}

# returns "unproductive" species of a reaction network
rn.unprods <- function(rn,mr=rn$mr,mp=rn$mp) {
  i <- which(colSums(mp)==0)  # unproductive reactions
  if (length(i)>0) {
    j <- which(rowSums(mr[,-i,drop=F])>0)  # directly productive species
    while (length(j)>0) {
      k <- which(colSums(mr[j,i,drop=F])>0)  # i[k] indirectly productive reactions
      if (length(k)==0) break
      j <- which(rowSums(mr[,i[k],drop=F])>0)  # indirectly productive species
      i <- i[-k]
    }
  }
  return(which(rowSums(mr[,i,drop=F])>0 | rowSums(mr)==0))
}

# returns a trimmed version of a reaction network rn (non reactive and consumed only species are recursively eliminated)
rn.trim <- function(rn) {
  mr <- rn$mr; mp <- rn$mp
  repeat {
    f <- T
    i <- rn.trans(NULL,mr,mp)
    if (length(i)>0) {
      f <- F
      j <- which(colSums(mr[i,,drop=F])>0)  # reactions using transient species
      mr <- mr[-i,-j,drop=F]; mp <- mp[-i,-j,drop=F]  # are eliminated as well as the species 
    }
    if (prod(dim(mr))==0) break
    i <- rn.unprods(NULL,mr,mp)
    if (length(i)>0) { f <- F; mr <- mr[-i,,drop=F]; mp <- mp[-i,,drop=F] }  # unused species are eliminated
    if (prod(dim(mr))==0) break
    i <- which(sapply(1:ncol(mr), function(k) all(mr[,k]==mp[,k])))
    if (length(i)>0) { f <- F; mr <- mr[,-i,drop=F]; mp <- mp[,-i,drop=F] }  # null reactions are eliminated
    if (prod(dim(mr))==0) break
    i <- which(duplicated(mr,MARGIN=2) & duplicated(mp,MARGIN=2))
    if (length(i)>0) { f <- F; mr <- mr[,-i,drop=F]; mp <- mp[,-i,drop=F] }  # duplicated reactions are eliminated
    if (prod(dim(mr))==0) break
    i <- which(duplicated(mr,MARGIN=1) & duplicated(mp,MARGIN=1))
    if (length(i)>0) { f <- F; mr <- mr[-i,,drop=F]; mp <- mp[-i,,drop=F] }  # complexes are reduced to the first species
    if (f || prod(dim(mr))==0) break
  }
  return(list(mr=mr,mp=mp))
}

# reassign the cardinalities of species in a reaction network rn
# range mass multiplier Rm sampled uniformly for each complex, mass conservation factor mcf
rn.recard <- function(rn,Rm=c(1,2),mcf=0) {
  for (i in 1:ncol(rn$mr)) {
    kr <- which(rn$mr[,i]>0); kp <- which(rn$mp[,i]>0)
    rn$mr[kr,i] <- 1; rn$mp[kp,i] <- 1
    Sr <- length(kr); Mr <- round(Sr*runif(1,Rm[1],Rm[2])) 
    Sp <- length(kp); Mp <- round(Sp*runif(1,Rm[1],Rm[2]))
    if (Mr*Mp!=0 && runif(1)<mcf)
      Mr <- Mp <- max(Sr,Sp,min(Mr,Mp)) # in case there is mass conservation
    if (Sr>0 & Mr>Sr) for (k in sample(Sr,Mr-Sr,replace=T)) rn$mr[kr[k],i] <- rn$mr[kr[k],i] + 1
    if (Sp>0 & Mp>Sp) for (k in sample(Sp,Mp-Sp,replace=T)) rn$mp[kp[k],i] <- rn$mp[kp[k],i] + 1
  }
  return(rn)
}

# returns the species (index) that belong to the closure of a set of species s (index) of a reaction network rn
# a is the set of active reactions (defaults to all reactions)
rn.closure <- function(rn,s,a=1:nrow(rn$mr)) {
  repeat {
    r <- which(colSums(rn$mr[,a,drop=F])==colSums(rn$mr[s,a,drop=F])) # triggered reactions
    ns <- which(rowSums(rn$mp[,a[r],drop=F])>0) # species produced by triggered reactions
    l <- length(s)
    s <- union(s,ns)
    if (length(s)==l) break # nothing new
  }
  return(sort(s))
}

# returns used species (consumed or produced by some triggered reaction) for a reaction network rn
# a set of species s (index) and a set of active reactions a (defaults to all reactions)
rn.used <- function(rn,s,a=1:nrow(rn$mr)) {
  r <- which(colSums(rn$mr[,a,drop=F])==colSums(rn$mr[s,a,drop=F])) # triggered reactions
  cs <- which(rowSums(rn$mr[,a[r],drop=F])>0) # species consumed by triggered reactions
  ps <- which(rowSums(rn$mp[,a[r],drop=F])>0) # species produced by triggered reactions
  return(sort(union(cs,ps)))
}

# returns a sub-network of a reaction network nr by closing the set of species s, unreactive species are trimmed
rn.close <- function(rn,s) {
  repeat {
    r <- which(colSums(rn$mr)==colSums(rn$mr[s,,drop=F])) # triggered reactions
    ns <- which(rowSums(rn$mp[,r,drop=F])>0) # species produced by triggered reactions
    l <- length(s)
    s <- union(s,ns)
    if (length(s)==l) break # nothing new
  }
  s <- sort(union(ns,which(rowSums(rn$mr[,r,drop=F])>0))) # support of triggered reactions
  srn <- list()
  if (!is.null(rn$w)) srn$w <- rn$w[s]
  srn$mr <- rn$mr[s,r]
  srn$mp <- rn$mp[s,r]
  return(srn)
}

# returns the basic sets of a reaction network rn (indexes to species s, reactions r, core reactions cr)
rn.basics <- function(rn) {
  b <- lapply(1:ncol(rn$mr), function(r) list(s=rn.closure(rn,rn.support(rn,r))))
  i <- 1
  while (i<length(b)) {
    j <- i + 1
    while (j<=length(b))
      if (setequal(b[[i]]$s,b[[j]]$s)) b[[j]] <- NULL else j <- j + 1
    i <- i + 1
  }
  for (i in 1:length(b)) {
    b[[i]]$cr <- b[[i]]$r <- unname(rn.supported(rn,b[[i]]$s))
    for (j in 1:length(b))
      if (i!=j && length(setdiff(b[[j]]$r,b[[i]]$r))==0) b[[i]]$cr <- setdiff(b[[i]]$cr,b[[j]]$r)
  }
  return(b)
}

# verifies that a reaction network is an organization (rn: reaction network, w: penalizing weights,
# inflow: given inflow (index of especies), destruct: penalize destruction (instead of creation)
# returns a set of needed extra inflow (empty set => organization)
rn.linp_org <- function(rn,w=NULL,inflow=NULL,destruct=F) {
  Ns <- nrow(rn$mr)  # number of species (rows)
  Nr <- ncol(rn$mr)  # number of reactions (columns)
  if(is.null(w)) {
    w <- rep(1,Ns) # penalization weights
    if (!is.null(inflow)) w[inflow] <- 0 # given inflow is not penalized
  }
  S <- cbind(diag(Ns),-diag(Ns),rn$mp-rn$mr) # stoichiometric matrix of rn with prepended creation and destruction reactions 
  f <- rep(0,Ns) # production of every species = 0, always possible because of additional creation and destruction reactions
  h <- rep(1.0,2*Ns+Nr); h[1:(2*Ns)] <- 0 # original reactions with rate>=1 (any positive), prepended reactions with rate>=0
  Cost <- rep(0,ncol(S)) # cost of processes is 0...
  Cost[1:Ns+Ns*destruct] <- w # except for creation (or destruction) reactions with cost w to be minimized
  rn.linp.r <<- linp(E=S,f,G=diag(2*Ns+Nr),h,Cost)
  # The unknown to be solved is the process column vector v = [vc,vd,vo] (creation/destruction/original network).
  # The result for v is stored in the variable rn.linp.r$X (vc is X[1:Ns], vd is X[Ns+(1:Ns)]).
  # The equations are S v = f (with f=0), G v >= h (i.e. [vc,vd,vo] >= [0,0,1] because G is the identity matrix).
  # Only the original network reactions are constrained to be strictly positive (here arbitrarily vo >= 1).
  # The cost is  Cost . vc  (dot product) or  Cost . vd  when destruction is penalized instead of creation.
  # The linear programming ideal result is every prepended creation reaction with rate 0 (minimal Cost 0).
  # It implies that there is at least one original process vector that can sustain every species (and even increase
  # them if destruction reactions have positive rate). In that case the reaction network is a proper organization.
  # If not it can be converted into one by adding the extra inflow of species with creation rate > 0.
  # Returns the indexes of species needed as extra inflow and overproduced species in the particular solution found.
  return(list(ifl=which(rn.linp.r$X[1:Ns]>0),ovp=which(rn.linp.r$X[Ns+1:Ns]>0)))
}

# verifies which species s (index set) of a reaction network rn (not necessarily an organization) are overproducible
# given an optional extra inflow (index set)
# returns a list of overproducible species (index set)
rn.overprod <- function(rn,s=1:nrow(rn$mr),inflow=integer(0)) {
  Ns <- nrow(rn$mr)  # number of species (rows)
  Nr <- ncol(rn$mr)  # number of reactions (columns)
  S <- cbind(-diag(Ns),diag(Ns)[,inflow],rn$mp-rn$mr) # stoichiometric matrix of rn with prepended reactions
  # of destruction of every species and optional creation of extra inflow species
  f <- rep(0,Ns) # flow vector constraint: production of every species = 0
  Cost <- rep(0,ncol(S)) # cost 0 for every reaction
  o <- rep(F,Ns) # overproducible status for species (initially False, until proven to be True)
  o[inflow] <- T # extra inflow species are known to be overproducible
  for (p in s) {
    if (o[p]) next  # already known overproducible species are skipped
    S[p,p] <- 1   # creation instead of destruction for p
    f[p] <- 1     # production of p = 1 instead of 0, if it is possible with creation rate 0, then it is overproducible
    Cost[p] <- 1  # cost for creation of p = 1 instead of 0
    r <- linp(E=S,f,Cost=Cost)
    # The unknown is the process vector v. The equations are S v = f with the inequality constraint v>=0 (rates can be 0).
    # Only the creation reaction for p is penalized in Cost (the rate should be 0 if p is overproducible).
    if (r$X[p]==0) {
      o[p] <- T  # no need of creation implies p is overproducible
      o[which(r$X[1:Ns]>0)] <- T # species with destruction rate > 0 are also overproducible
    }
    S[p,p] <- -1; f[p] <- 0; Cost[p] <- 0  # original destruction reaction and zero values for next iteration
  }
  return(which(o))
}

# production connections between species within a reaction network rn
rn.pconect <- function(rn) {
  r <- rn$mr>rn$mp
  p <- rn$mp>rn$mr

}

rn.example <- function() {
  rn <<- rn.testrn(1)
  cat("*** farm:\n")
  rn.display(rn)
  cat("overproducible species:",rownames(rn$mr)[rn.overprod(rn)],"\n")
  ss <- sapply(c("grass","cows","milk","worms","dung"), function(s) match(T,rownames(rn$mr)==s))  # the index of named species
  srn <<- rn.sub(rn,s=setdiff(1:nrow(rn$mr),ss))
  cat("*** dairyless farm:\n")
  rn.display(srn)
  l <- rn.linp_org(srn)$ifl
  cat("needed inflow:",rownames(srn$mr)[l],"\n")
  cat("overproducible species:",rownames(srn$mr)[rn.overprod(srn)],"\n")
}

rn.test <- function(i=3,rn=rn.testrn(i)) {
  rn <<- rn
  rn.display(rn)
  r <<- rn.linp_org(rn)
  cat("needed inflow:",rn$sid[r$ifl],"\n")
  cat("overproducible species:",rn$sid[rn.overprod(rn)],"\n")
}
