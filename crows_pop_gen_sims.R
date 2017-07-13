###
# DO SOME POPULATION GENETIC SIMULATIONS
###
options(scipen = 6)
# REFERENCES:
#(1) Charlesworth and Charlesworth 2008 Elements of Evolutionary Genetics
#(2) Smeds et al., 2016 Genome Research 26:1211-1218
#(3) Ellegren 2013 Annual Review of Ecology, Evolution and Systematics 44:239-259
#(4) Nam et al., 2010 Genome Biology 11: R68
args<-commandArgs(trailingOnly = TRUE)
# Sample XX "SNPs" to make up an SFS
nsnps <- as.numeric(args[1])
# How many times to run each simulation
nsims <- as.numeric(args[2])
cat("N SNPs = ", nsnps, " - ","N Sims = ", nsims,"\n")

source("~/RData/RScripts/sfs_stats.R")
se<-function(x){sqrt(var(x,na.rm=TRUE))/sqrt(length(x))}

# RUN SIMULATIONS
# Simulate Tajima's D and pi for different sample sizes and Ne.
# Assume a neutral population at mutation-drift balance. 
# Take different effective population sizes
Ne<-c(1000,10000,20000,100000,200000,1000000,2000000)
# Mutation rate (same as for humans)
MU<-c((2.3*10^-9),(2.21*10^-9),(1.91*10^-9),(1.23*10^-9))# REFs: (2), (4)
# Take different numbers of diploid individuals from the population
N<-c(5,10,50,100)

# Make vectors to store data
Nes<-vector(length=nsims*length(Ne)*length(N)*length(MU))
Ns<-vector(length=nsims*length(Ne)*length(N)*length(MU))
mus<-vector(length=nsims*length(Ne)*length(N)*length(MU))
tajDs<-vector(length=nsims*length(Ne)*length(N)*length(MU))
pis<-vector(length=nsims*length(Ne)*length(N)*length(MU))
alphs<-vector(length=nsims*length(Ne)*length(N)*length(MU))
i <- 1
for(ne in Ne){
  for(mu in MU){
    for(n in N){
      for(sim in 1:nsims){
        cat(paste("Ne = ",ne,sep=""),"-",paste("mu = ",mu,sep=""),"-",
            paste("N = ",n,sep=""),"-",paste("sim = ",sim,sep=""),"\n")
        SFS<-vector(length=(2*n)+1)
        for(snp in 1:nsnps){
          # Sample an allele frequency from the betabinomial distribution.
          # Frequency of A allele
          alpha<-4*ne*mu
          pA<-rbeta(n=1,shape1=4*ne*mu,shape2=4*ne*mu) #REF: (1)
          # Ge the number of genotypes for sample of n individuals assuming HWE
          # 1 = p^2 + 2pq + q^2
          gen_typs<-rmultinom(n=1,size=n,prob=c((pA^2),(2*pA*(1-pA)),((1-pA)^2)))
          # Get the number of A alleles in the sample of n diploid individuals
          # This is 2*Homozygotes(A) + Heterozygotes(Aa)
          ac<-sum(2*gen_typs[1],gen_typs[2])
          # SFS index is then ac + 1
          SFSi<-ac+1
          SFS[SFSi]<-SFS[SFSi]+1
        }
        Nes[i]<-paste("Ne = ",ne,sep="")
        Ns[i]<-paste("N = ",n,sep="")
        tajDs[i]<-tajD_sfs(SFS)
        pis[i]<-pi_sfs(SFS)
        alphs[i]<-alpha
        mus[i]<-mu
        i<-i+1
      }
    }
  }
}

sim_dat<-data.frame(Ne=Nes,N=Ns,tajd=tajDs,pi=pis,alpha=alphs,mu=mus)
# Save the simulation data:
write.table(sim_dat,
            paste("crows_pop_gen_sims_ns",nsnps,"_s",nsims,".tab",sep=""),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")