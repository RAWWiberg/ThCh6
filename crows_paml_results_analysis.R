###
# Analysis of the crow codeml results
#
#
###
# Clear all
rm(list = ls(all=TRUE))

#
# Load libraries #####
library(ggplot2)
library(stringr)
source("~/RData/RScripts/ggplot_theme.R")
library(reshape2)
library(dplyr)
library(plyr)
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library(qvalue)
# ####

#
species8<-c("3sp_Ccornix","Ccorx","Chaw","Cfru","Csple","Ctas","Cmon","Cdau","Tgut")
species7<-c("3sp_Ccornix","Ccorx","Cfru","Csple","Ctas","Cmon","Cdau","Tgut")
species5<-c("3sp_Ccornix","Cfru","Csple","Cmon","Cdau","Tgut")
# Load data ####
nspecies<-"8species"
seqs<-"_conv"
A_data<- read.table(paste("crows_paml/",
                          nspecies,
                          "_A_results_pamlinfiles",seqs,".txt",sep=""),
                    sep=",",header = FALSE)
A_data$V1<-gsub(" ","",A_data$V1)
head(A_data)
# Add header to dataset
if(nspecies=="8species"){
  # 8 species dataset: 3sp_Ccornix,Ccorx,Chaw,Cfru,Csple,Ctas,Cmon,Cdau,Tgut ###
  colnames(A_data)<-c("gene","stop","N","S",
                      paste(species8,"_dN",sep=""),
                      paste(species8,"_dS",sep=""),
                      paste(species8,"_w",sep=""),
                      "lnL","np","model")
  # Load the 8species complete genes
  complg<-read.table(paste("crows_paml/",
                   nspecies,
                   "_complete_genes_pamlinfiles.list",sep=""),
             sep=",",header = FALSE)
  colnames(complg)<-c("gene")
  # Subset
  A_data<-A_data[A_data$gene %in% complg$gene,]
}else if(nspecies=="7species"){
  # 7 species dataset: 3sp_Ccornix,Cdau,Cfru,Csple,Cmon,Ctas,CcorxTgut ###
  colnames(A_data)<-c("gene","stop","N","S",
                      paste(species7,"_dN",sep=""),
                      paste(species7,"_dS",sep=""),
                      paste(species7,"_w",sep=""),
                      "lnL","np","model")
}else if(nspecies=="5species"){
  # 5 species dataset: 3sp_Ccornix,Cdau,Cfru,Csple,Cmon,Tgut ###
  colnames(A_data)<-c("gene","stop","N","S",
                      paste(species5,"_dN",sep=""),
                      paste(species5,"_dS",sep=""),
                      paste(species5,"_w",sep=""),
                      "lnL","np","model")
}

# ####
Null_data<- read.table(paste("crows_paml/",
                             nspecies,
                             "_Null_results_pamlinfiles",seqs,".txt",sep=""),
                       sep=",",header = FALSE)
Null_data$V1<-gsub(" ","",Null_data$V1)
# Add header to dataset
if(nspecies=="5species"){
  # 5species dataset:
  colnames(Null_data)<-c("gene","stop","N","S","dN",
                         "dS","w","lnL","np","model",
                         paste(species5,"_dN",sep=""),
                         paste(species5,"_dS",sep=""))
  }else if(nspecies == "7species"){
    # 7species dataset:
    colnames(Null_data)<-c("gene","stop","N","S","dN",
                           "dS","w","lnL","np","model",
                           paste(species7,"_dN",sep=""),
                           paste(species7,"_dS",sep=""))
  }else if(nspecies == "8species"){
    # 8species dataset:
    colnames(Null_data)<-c("gene","stop","N","S","dN",
                           "dS","w","lnL","np","model",
                           paste(species8,"_dN",sep=""),
                           paste(species8,"_dS",sep=""))
    # Subset
    Null_data<-Null_data[Null_data$gene %in% complg$gene,]
  }

# ####
nrow(A_data)
nrow(Null_data)

head(A_data)
head(Null_data)
#
# Perform LRT ####
A_data$p <- vector(length=nrow(A_data))
for (g in A_data$gene){
  np1 <- A_data$np[A_data$gene == g]
  np2 <- Null_data$np[Null_data$gene == g]
  lnL1 <- A_data$lnL[A_data$gene == g]
  lnL2 <- Null_data$lnL[Null_data$gene == g]
  LR <- 2*(lnL1 - lnL2)
  df <- np1-np2
  A_data$p[A_data$gene == g]<-pchisq(LR,df,lower.tail = FALSE)
  Null_data$p[Null_data$gene == g]<-pchisq(LR,df,lower.tail = FALSE)
}
# ####

#
# Subset from A_data those with significant LRT ####
# Write p-values to a list
write.table(A_data$p,paste("crows_paml/",nspecies,seqs,"_sep_A_data_pval.list",sep=""),
            row.names = FALSE,quote = FALSE,col.names = FALSE)
# Load p-values to make q-values
pvals<-scan(paste("crows_paml/",nspecies,seqs,"_sep_A_data_pval.list",sep=""))
qobj <- qvalue(pvals)
plot(qobj)
A_data$qvalues <- qobj$qvalues
Null_data$qvalues <- qobj$qvalues
head(A_data)

# Which are the "significant" genes (q-value < 0.05)
sig_A_data <- A_data[A_data$qvalues < 0.05,]
# Remove the spaces in the gene names
sig_A_data$gene <- gsub(" ","",sig_A_data$gene)
# How many
nrow(sig_A_data)
sig_A_data[,1]
# Show the data
head(sig_A_data)
tail(sig_A_data)
# nrow(sig_A_data[sig_A_data$p < bonf,])
sig_A_data
A_data[A_data$gene=="DTNBP1 ",]

beak_genes<-c("POU1F1","IGF2R","BMP4","")

head(A_data[A_data$p < 0.05,])
sig_A_data<-A_data[A_data$p < 0.05,]
sig_A_data[which(sig_A_data$gene %in% genes),]

# Write significant results to a table
write.csv(sig_A_data,
          paste("~/PhD/crows_project/results/paml/",
                nspecies,
                seqs,
                "_sig_A_data.tab",
                sep=""),
          row.names = FALSE,quote = FALSE)

# ####








