# DEFINE FUNCTIONS FOR GETTING STATS FROM SFS

# Source Daniel Halligan's SFS functions
#source("RScripts/polymorph.R")

# References:
# (1) Tajima 1989 Genetics 123:585-595
# (2) Charlesworth & Charlesworth 2008 Elements of Evolutionary Genetics
# (3) Daniel L. Halligan pers. comm.
# (4) Nielsen et al., 2012 PLoS ONE 7:e37558
# (5) Korneliussen et al., 2013 BMC Bioinformatics

# Define a "testing" SFS
SFS<-c(4953,19,11,10,12,5,10,9,11,28,4932)

# Function to "Fold" the SFS
# Will sum the first half of the SFS to it's reverse half
# If the SFS is an odd length the middle bin is halved.
# The folded/end half of the SFS is set to 0.
foldSFS<-function(sfs){
  l<-length(sfs)
  c1<-floor(l/2)
  c2<-ceiling(l/2)
  fSFS<-sfs+rev(sfs)
  # If c1 and c2 are not the same the length of the SFS is odd
  # So halve the middle bin 
  if(c1!=c2){sfs[c2]<-sfs[c2]/2}
  fSFS[(c2+1):length(fSFS)]<-0
  return(fSFS)
}


# a1 = 1/1 + 1/2 + 1/3 + 1/4 ... 1/(n-1), where 
# n = the number of sequences compared (2*N for diploids) 
a1_sum<-function(n){return(sum(1/(1:(n-1))))} # REF: (1)
# a2 = 1/1^2 + 1/2^2 + 1/3^2 + 1/4^2 ... 1/((n-1)^2)
a2_sum<-function(n){return(sum(1/(1:(n-1))^2))} # REF: (1)

###
# Get pi from SFS
# This is possible if we assume that each site is biallelic
# REF: (1) (5)
# Called k_hat in (1)
# pi = w = 4Neu (under the neutral model)
###
# k is sum(h); where
# h = (n*(1-sum(x^2)))/(n-1)
# where x is the frequency of nucleotides
pi_sfs<-function(sfs,persite=TRUE){
  # The number of aligned sequences
  n<-length(sfs)-1
  if(n<=1){return(NA)}
  i<-(n:1)
  pi<-sum(i*(n-i)*sfs[1:n])/choose(n,2)
  if(persite){
    return(pi/sum(sfs))
  }else{
    return(pi)
  }
}
pi_sfs(SFS)
pi_sfs(foldSFS(SFS))
###
# Get S (Watterson's theta) from SFS
# REF: (1)
# Called M_hat in (1) or S/a1
###
w_sfs<-function(sfs,persite=TRUE){
  n<-length(sfs)-1
  if(n<=1){return(NA)}
  # How many polymorphic sites
  S<-sum(sfs[2:n])
  a1<-a1_sum(n)
  w<-S/a1
  if(persite){
    return(w/sum(sfs))
  }else{
    return(w)
  }
}
pi_sfs(SFS)
pi_sfs(foldSFS(SFS))

###
# Get tajima's D from SFS
# REF: (1)
###
# TajD = (pi - (S/a1))/sqrt(V); where
# V = var(pi - (S/a1)) # REF: (1)
# V is given by:
# (e1*S)+(e2*S*S-1); where
# e1 = c1 - (1/a1) # REF: (1)
# and
# e2 = c2 / (a1^2 + a2) # REF: (1)
# where;
# b1 = n+1/(3*(n-1)) # REF: (1)
# b2 = (2*((n^2)+n+3))/(9*n*(n-1)) # REF: (1)
# c1 = b1 - (1/a1) # REF: (1)
# c2 = b2 - ((n+2)/a1*n)+(a2/a1^2) # REF: (1)

tajD_V<-function(n,S){
  a1<-a1_sum(n)
  a2<-a2_sum(n)
  
  b1<-(n+1) / (3*(n-1))
  b2<-(2*(n^2 + n + 3)) / (9*n*(n-1))
  
  c1<-b1 - (1/a1)
  c2<-b2-((n+2)/(a1*n))+(a2/a1^2)
  
  e1<-c1/a1
  e2<-c2/((a1^2)+a2)
  V<-(e1*S)+(e2*S*(S-1))
  return(V)
}
tajD_sfs<-function(sfs){
  n<-length(sfs)-1
  S<-sum(sfs[2:n])
  D<-(pi_sfs(sfs,persite=FALSE)-w_sfs(sfs,persite=FALSE))/sqrt(tajD_V(n,S))
  return(D)
}
