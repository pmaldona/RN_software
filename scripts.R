
##########################################
# Example scripts for using the software #
##########################################

# packages needed: limSolve, ggplot2, reshape2, cowplot, magrittr, viridisLite, jsonlite

# R files to be loaded:
# rg.R (random generation of reaction networks)
# rn.R (structural analysis and manipulation of reaction networks including linear programming)
# sm.R (dynamical simulation of reaction networks)
# pert.R (state and flow perturbations)
# ev.R (evolution potential of reaction networks calculated using many dynamical simulations of random perturbations)

library("jsonlite")
library("parallel")
source("ev.R")  # all other R files are loaded by ev.R

#-----------------------------------------------------------------------------
# A script for generating random reaction networks and testing their dynamics
# Use: scr.test()
# It uses scritpts defined below
#

scr.test <- function() {
  scr.gen() # a global variable rn with a random generated reaction network will be created
  scr.simul(rn) # will simulate the dynamics of rn and display the end result and time graphs
}

#-------------------------------------------------
# Scripts for generating random reaction networks
#

scr.gen <- function(case=1) {
  rn <<- switch(case,
                rg.g1(Nr=20,Ns=10,extra=.5), # case 1: Nr reactions, Ns species
                rg.bg(Nr=20,Rs=1:3,Rm=c(1,2),dsf=.2) # case 2: Nr reactions, number of species depends on Rs and dsf
  )
  rn <<- rn.merge(rn) # null or redundant reactions are filtered out, we may end with less reactions and species
  trn <<- rn.trim(rn)  # a more elaborate filtering that eliminates also non reactive and obvious transient species
}

# a general script with a more complex generator
scr.genrn2 <- function(Nr=100,Ns=Nr) {
  
  # (1) reaction patterns to be used:
  mP <- rbind(c(0,0,0,1,1,50),c(0,0,0,1,2,20),c(0,0,0,2,1,20),c(0,0,0,2,2,10))
  # each pattern is a 6-tuple of species numbers in different roles and weight of the pattern
  # the roles are catalizer, mixed-consumed, mixed-produced, pure-consumed, pure-consumed
  # c(0,0,0,1,1,50) means a reaction with one pure-consumed and one pure-produced species, like a => b
  # c(0,0,0,1,2,50) means a reaction with one pure-consumed and 2 pure-produced species, like a => b + c
  # c(1,1,0,0,0,33) would mean a reaction with one catalyzer and one mixed-consumed species, like a + 2b => a + b
  # each species belonging to a pattern will be a different one when the reaction is generated
  
  # (2) number of reactions to be generated: Nr, a parameter defined by default to be 100
  
  # (3) species (names) to be used:
  S <- paste0("s",1:Ns)
  # S = ["s1","s2","s3",... , "s100"]
  # By default we are using the same number of species as reactions (parameter Ns). Names can be any string.
  # Species are selected at random, many species may not appear in the resulting reaction network.
  
  # (4) cardinalities of species:
  C <- 1:3
  # C = [1,2,3], here these three values are the only ones available for the generation of reactions
  # so s1 => 3s2 (cardinalities 1 an 3) can be generated, but not 5s1 => 4s2 (cardinalities 5 and 4).
  # If there is only one cardinality, mixed species are generated with a second cardinality adding 1.
  # Non integer positive cardinalities can also be used.
  
  # (5) pattern generation function:
  P <- function() rg.rsrt(m=mP)
  # the function uses the patterns just defined, it chooses a pattern at random according to weight
  # we need this function to be given as a parameter to the the generating function rg.rn
  
  # (6) distribution of species (D$s) and cardinalities (D$c)
  D <- list(s=length(S)/(1:length(S))+length(S)/5,c=length(C)/(1:length(C)))
  # D$s should be a matrix of dimensions Ns X 5 were each of the five columns represents the weight of the Ns species
  # according to the corresponding role in the same order as defined in the reaction patterns (one column per role).
  # If D$s has fewer than 5 columns, it is completed by repeating the last column. Here we have only a vector,
  # it is interpreted as a column an the matrix is completed as described (all the columns are the same).
  # D$c is a vector with the weigth of each cardinality.
  # Here both D$s and D$c give much more weight to first values than last, s1 will appear much more frequently than
  # s100, and cardinality 1 more than cardinality 3. For D$s it is important to define highly frequent species so that
  # the generated reaction network form a connected unit, specially if there are not many reactions.
  
  # (7) generation of a random network using the defined parameters:
  rn <<- rg.rn(Nr=Nr,S=S,C=C,P=P,D=D)
  # the result rn is created as a global variable that will persist after the end of this script.
  # There is an additional parameter that can be used, rep=T, which indicates that we accept the repeated generation
  # of a reaction (by default rep=F).
  
  # (8) trimming the generated reaction network:
  trn <<- rn.trim(rn)
  # the result trn is created as a global variable that will persist after the end of this script.
  # non reactive species are eliminated as well as obviously transient species (by statical analysis) and the reactions
  # that use them, repeated and null reactions are also eliminated.
  
  return(trn)
}

# a modification of the last script to generate some reactions with exclusive catalyzers
scr.genrn3 <- function(Nr=100,Ns=Nr) {
  
  # (1) reaction patterns to be used:
  mP <- rbind(c(1,0,0,1,1,10),c(0,0,0,1,1,40),c(0,0,0,1,2,20),c(0,0,0,2,1,20),c(0,0,0,2,2,10))
  # we add a first pattern with a catalyzer with weight 10, the former first pattern weight is reduced from 50 to 40
  
  # (3) species (names) to be used:
  S <- c(paste0("c",1:2),paste0("s",1:Ns))
  # S = ["c1","c2","s1","s2","s3",... , "s100"], we added species "c1" and "c2" to be used as catalyzers.
  
  # (4) cardinalities of species:
  C <- 1:3
  
  # (5) pattern generation function:
  P <- function() rg.rsrt(m=mP)
  
  # (6) distribution of species (D$s) and cardinalities (D$c)
  D <- list()
  D$s <- cbind(c(1,1,rep(0,Ns)), c(0,0,(length(S)-2)/(1:(length(S)-2))+(length(S)-2)/5))
  D$c <- length(C)/(1:length(C))
  # D$s is now defined as a matrix of 2 columns (that will be completed repeating the last column 3 times).
  # The first column describes the species that will be chosen to be catalizers, only the first two (c1 and c2)
  # with the same weight. The second column stands for all the others choices, which is every species but the
  # first two (c1 and c2 will never play a different role than catalyzers).
  
  # (7) generation of a random network using the defined parameters:
  rn <<- rg.rn(Nr=Nr,S=S,C=C,P=P,D=D)
  
  # (8) trimming the generated reaction network:
  trn <<- rn.trim(rn)
  
  return(trn)
}

# a script to display generated reaction networks (by default it displays the global variable trn)
scr.rndisp <- function(rn=trn) {
  
  # (1) displaying the number of species and reactions
  cat("number of species:",nrow(rn$mr),"; number of reactions:",ncol(rn$mr),"\n")
  # The names of the species and reactions is rownames(rn$mr) and colnames(rn$mr).
  # The reaction network is represented as two matrices, one for reactants (rn$mr) and the other for products (rn$mp).
  
  # (2) displaying the reactions in a readable format
  rn.display(rn)
  # The reactions will be displayed in the console in Antimony like format.
  # If the parameter file="name_of_file" is defined, then a file is created.
  # Files can be read using rn.read("name_of_file").
  
  # (3) using linear programming to assess the reaction network
  r <- rn.linp_org(rn)
  cat("needed inflow:",rn$sid[r$ifl],"\n")
  cat("overproducible species:",rn$sid[rn.overprod(rn)],"\n")
  # If the needed inflow is void, then n is an organization (but maybe a dynamically unstable one).
  
}

#-------------------------------------------------------------
# Script to test the dynamic simulation of a reaction network
#

scr.simul <- function(rn=sm.genrn(12),n=1000,dt=.2,e=.2) {
  
  # (1) by default we are using a small random generated reaction network produced by a different algorithm
  # we display the network
  cat(nrow(rn$mr),"species,",ncol(rn$mr),"reactions:\n")
  rn.display(rn)
  
  # (2) using linear programming to statically assess the reaction network
  o <- rn.linp_org(rn)
  cat("needed inflow {",rownames(rn$mr)[o$ifl],"}\n")
  cat("overproduced if inflow {",rownames(rn$mr)[o$ovp],"}\n")
  cat("overproducible {",rownames(rn$mr)[rn.overprod(rn)],"}\n")
  
  # (3) a mass action kinetics simulation for n iterations, dt time step and e existence threshold 
  sm <<- sm.maksim(rn,n=n,dt=dt,e=e)
  
  # (4) we display the final state (the closed set abstraction of the actual final state)
  s <- sm.final(rn,sm)
  cat("final closed abstract state {",rownames(rn$mr)[s],"}\n")
  rs <- rn.supported(rn,s) # reactions supported by the final state
  cat("final closed abstract activity: {",colnames(rn$mr)[rs],"}\n")
  cat("missing reactions: {",colnames(rn$mr)[(1:ncol(rn$mr))[-rs]],"}\n")
  
  # (5) graphic displaying of the simulation
  gc()  # garbage collection... all unused memory space in R is recovered
  g <- sm.display(sm,simple=F) # the graphs are generated (simple=T would generate a simpler graph)
  print(g) # the graphs are displayed
  Sys.sleep(1) # to give time so the graphs are displayed before exiting the current script
}

#--------------------------------------------------------------------------
# Scripts related to testing the evolutive potential of a reaction network
#

# testing of a reaction network using systematically the void as starting point and perturbations of any size
scr.evol <- function(rn=scr.genrn(),M=5000,n=1000) {
  
  # (1) calculates the evolutive potential of rn (by default a random generated reaction network)
  # for M random perturbations and n iterations of a mass action kinetics to reach the final state
  e <<- ev.evol(rn=rn,M=M,n=n)
  
  # (2) analyses the raw results of the first step
  ae <<- ev.analyse(e)
  
  # (3) displays some distributions and indicators over the analysis
  ev.disp(ae)
  
}

# testing of a reaction network using random starting points and perturbations of size 1
# the script is designed to complete the results of a previous test stored in global variable e
# if e is not available a new random reaction network is used instead
scr.evol2 <- function(e=globalenv()$e, rn = if (!is.null(e)) e$rn else scr.genrm(),M=5000,n=1000) {
  
  # (1) we define a distribution for the size of perturbations, in this case all perturbations will be of size 1
  P <- function(n) c(1,rep(0,n-1))
  # the size of a perturbation is the number of *reactions* we add
  
  # (2) we add M more perturbations to e starting from randomly selected final states previously generated
  e2 <<- ev.evol(e,rn,P=P,M=M,n=n,systematic=F)
  # by default systematic=T meaning that we are always perturbing the same state, by default the closure of the void
  # systematic=F means a random selection of the state to be perturbed. The states are the ones that were found as
  # final states of previous perturbations. New perturbation are incrementally added to the already calculated ones
  # in e. The parameter P is used to indicate the distribution of perturbation size.
  
}

# working script to generate a random reaction network and perturbation/simulation random walks
# if e is provided new random walks or steps are added to e
# if rn is provided that reaction network is used instead of generating a random one
# the random walks to be created or completed are defined by range w, by default 1:10
# by default the number of perturbation and simulation steps l is set to 10
# cutoff is the concentration threshold for a species to be reactive and n is the number of simulation steps 
# the result is available in global variable e (the value is also returned by the function)
scr.gen_and_pert <- function(e=NULL,rn=NULL,w=1:10,l=10,cutoff=.1,n=5000) {
  if (is.null(e)) {
    if (is.null(rn)) {
      rn <- rg.g1(Nr=100,Ns=100,extra=.5) # random reaction network with Nr reactions, Ns species, .5 extra species
      rn <- rn.merge(rn) # null or redundant reactions are filtered out
    }
    e <- pert.start(rn) # the starting structure to store the random walks to be generated
  }
  for (i in w) { # for each random walk
    if (i>length(e$rw)) # this is a new random walk
      e$rw[[i]] <- list(t=NULL,f=NULL,s=NULL,p=NULL,c=NULL,a=NULL,u=NULL) # matrices are created to store the steps in columns
    if (is.null(e$rw[[i]]$f)) { # this is a void random walk (0 steps)
      s <- e$rn$mr[,1]*0 # the current state is zeroed
      f <- pert.randomize(e$rn$mr[1,]*0 + 1) # the flow vector is randomized (each random walk has a different f)
    }
    else { # this is a random walk with a number of steps already accumulated
      j <- ncol(e$rw[[i]]$c)  # the number of steps up to now
      s <- e$rw[[i]]$c[,j] # the exploration continues from the last convergence state in the random walk
      f <- e$rw[[i]]$f[,j] # the last flow vector is conserved
    }
    for (j in 1:l) { # for each j step in random walk i
      s <- s*(s>cutoff) # species of the current state under the cutoff threshold are zeroed
      # if (all(s>0)) break # no more species to add... the random walk has finished before the l step
      e$rw[[i]]$s <- cbind(e$rw[[i]]$s,unname(s)) # the current state is stored in the random walk
      # start perturbations:
      s <- pert.delta(s,d=1,nmin=1,sigma=.5) # a delta perturbation is applied to current state
      # end perturbations, start simulation:
      st <- system.time(
        cs <- pert.simul(e$rn,s,f,cutoff=cutoff,n=n), # the perturbed state is simulated reaching a convergence state
        gcFirst = FALSE
      )
      # end simulation
      e$rw[[i]]$t <- c(e$rw[[i]]$t,st[1]) # the time elapsed in the dynamic simulation is stored in the random walk
      e$rw[[i]]$f <- cbind(e$rw[[i]]$f,unname(f)) # the flow vector is stored in the random walk
      e$rw[[i]]$p <- cbind(e$rw[[i]]$p,unname(s)) # the perturbed current state is stored in the random walk
      if (is.null(cs)) cs <- s # if something went wrong with the simulation, we just keep the initial current estate
      e$rw[[i]]$c <- cbind(e$rw[[i]]$c,unname(cs)) # the convergent state is stored in the random walk
      a <- pert.abstract(e$rn,cs,f,cutoff=cutoff) # boolean abstraction of the convergent state (closure)
      e$rw[[i]]$a <- cbind(e$rw[[i]]$a,unname(a+0)) # the abstraction is stored
      e$rw[[i]]$u <- cbind(e$rw[[i]]$u,pert.used(e$rn,a,f)+0) # a second abstraction is stored (used species)
      s <- cs # the new current state is the convergent state
    }
  }
  e <<- e
}

# saves variable e in json format into a file
scr.save <- function(file="rw.json") {
  write(toJSON(e,digits=NA),file)
}

# loads a json file and deparses it into variable e and returns it
scr.load <- function(file="rw.json") {
  e <- fromJSON(file,simplifyDataFrame=F)
  rownames(e$rn$mr) <- rownames(e$rn$mp) <- e$species
  colnames(e$rn$mr) <- colnames(e$rn$mp) <- e$reactions
  e <<- e
}

# generates n random walks and stores each random walk in a file
scr.random_walk <- function(n=10,file="rndw") {
  if (n>0) for (i in 1:n) {
    rn <- sm.genrn(12) # a new random network
    scr.gen_and_pert(rn=rn,w=1:10,l=10,cutoff=.1,n=5000) # a random walk for rn
    scr.save(file=paste0(file,i,".json")) # result stored in files rndw1.json, rndw2.json, etc. according to index i
  }
}

# generates a dataframe of combinations of parameters
# reaction networks: nr number of reactions, ns number of species, w random walks of l steps, N different tries
# conditions: n simulation iterations with the same reaction networks
scr.pcomb <- function(rn=cbind(nr=c(50,100,250),ns=c(50,100,250),w=100,l=100),N=1:5,n=c(1000,2500,5000)) {
  d <- data.frame(rn)
  d <- merge(d,data.frame(N))
  d$id <- 1:nrow(d)
  d <- merge(d,data.frame(cbind(n=n,cond=1:length(n))))
  return(d)
}

# generates a list of reactions networks according to parameters in dataframe p
# if the same identity is used several times, the first occurrence in p is the defining one
scr.pcomb.genr <- function(p,rng=1) {
  rngen <- function(nr,ns,rng) {  # generates a random network using generator rng
    if (rng==1) {
      patterns <- rbind( c(1,0,0,0,1,10),  # creation
                         c(1,0,0,1,1,80),  # transformation
                         c(1,0,0,1,0,10) ) # destruction
      rn <- rg.rn(Nr=nr,S=paste0("s",1:ns),C=1, P = function() rg.rsrt(patterns))
    }
    else if (rng==2) {
      rn <- rg.g1(Nr=nr,Ns=ns,extra=.2, dist=function(x) -((Ns-1)/2*x/5)^2, pr=log(100), pp=log(100))
      rn <- rn.merge(rn) # null or redundant reactions are filtered out
    }
    else if (rng==3) {
      rn <- sm.genrn(nr)
    }
  }
  u <- unique(p$id) # the unique identities of reaction networks to be used
  rnl <- lapply(match(u,p$id), function(i) rngen(p$nr[i],p$ns[i],rng))
  names(rnl) <- u
  return(rnl)
}

# batch function to create random walks according to parameters in dataframe p
scr.batch <- function(p=scr.pcomb(),rnl=scr.pcomb.genr(p,rng=1),name="batch") {
  d <- cbind(id=1:length(rnl),t(sapply(rnl,function(e) c(gnr=ncol(e$mr),gns=nrow(e$mr)))))
  p <- merge(p,as.data.frame(d),by="id")
  for (i in 1:nrow(p)) {
    st <- system.time(scr.gen_and_pert(rn=rnl[[p$id[i]]],w=1:p$w[i],l=p$l[i],cutoff=.1,n=p$n[i]))
    d <- data.frame(id=p$id[i],cond=p$cond[i],time=st[1])
    write.csv(d,sprintf("%s%04i_%02i.csv",name,p$id[i],p$cond[i]),row.names=F)
    scr.save(file=sprintf("%s%04i_%02i.json",name,p$id[i],p$cond[i]))
  }
  p$time <- 0
  for (i in 1:nrow(p)) {
    d <- read.csv(sprintf("%s%04i_%02i.csv",name,p$id[i],p$cond[i]))
    p$time[i] <- d$time
  }
  write.csv(p,sprintf("%s.csv",name),row.names=F)
}

# batch function to create random walks according to parameters in dataframe p
scr.batch.parallel <- function(p=scr.pcomb(),rnl=scr.pcomb.genr(p,rng=1),name="batch") {
  d <- cbind(id=1:length(rnl),t(sapply(rnl,function(e) c(gnr=ncol(e$mr),gns=nrow(e$mr)))))
  p <- merge(p,as.data.frame(d),by="id")
  p_list <- split(p,1:nrow(p))
  
  int.batch <- function(param,name,rnl){
    st <- system.time(scr.gen_and_pert(rn=rnl[[param$id]],w=1:param$w,l=param$l,cutoff=.1,n=param$n))
    d <- data.frame(id=param$id,cond=param$cond,time=st[1])
    print(d)
    write.csv(d,sprintf("%s%04i_%02i.csv",name,param$id,param$cond),row.names=F)
    scr.save(file=sprintf("%s%04i_%02i.json",name,param$id,param$cond))
  }
  mclapply(p_list,FUN = function(x){int.batch(x,name,rnl)},mc.cores = 24, mc.preschedule = T)
  p$time <- 0
  for (i in 1:nrow(p)) {
    d <- read.csv(sprintf("%s%04i_%02i.csv",name,p$id[i],p$cond[i]))
    p$time[i] <- d$time
  }
  write.csv(p,sprintf("%s.csv",name),row.names=F)
}
