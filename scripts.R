
##########################################
# Example scripts for using the software #
##########################################

# packages needed: limSolve, ggplot2, reshape2, cowplot, magrittr, viridisLite

# R files to be loaded:
# rg.R (random generation of reaction networks)
# rn.R (structural analysis and manipulation of reaction networks including linear programming)
# sm.R (dynamical simulation of reaction networks)
# ev.R (evolution potential of reaction networks calculated using many dynamical simulations of random perturbations)

source("ev.R")  # all other R files are loaded by ev.R

#-------------------------------------------------
# Scripts for generating random reaction networks
#

# a general script
scr.genrn <- function(Nr=100,Ns=Nr) {

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

# a modification of the first script to generate some reactions with exclusive catalyzers
scr.genrn2 <- function(Nr=100,Ns=Nr) {
  
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

#--------------------------------------------------------------------------
# Scripts related to testing the evolutive potential of a reaction network
#

scr.evol <- function(rn=scr.genrn(),M=5000,n=1000) {

  # (1) calculates the evolutive potential of rn (by default a random generated reaction network)
  # for M random perturbations and n iterations of a mass action kinetics to reach the end state
  e <<- ev.evol(rn=rn,M=M,n=n)
  
  # (2) analyses the raw results of the first step
  ae <<- ev.analyse(e)
  
  # (3) displays some distributions and indicators over the analysis
  ev.disp(ae)
  
}
