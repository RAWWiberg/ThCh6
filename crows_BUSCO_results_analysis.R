###
# Assessing the crow BUSCO results
# Last Modified: 29.04.2017
###
#clear all
rm(list = ls(all=TRUE))

# Load libraries #####
library(ggplot2)
library(stringr)
source("~/RData/RScripts/ggplot_theme.R")
library(reshape2)
library(dplyr)
library(plyr)
# ####
res<-data.frame()
species<-c("3sp_Ccornix","Ccorx","Cdau","Cfru",
           "Chaw","Cmon","Csple","Ctas","Cbrach",
           "Ckub","Cwoo","Cmone","Tgut")
# read in the data
for(sp in species){
  dat<-read.table(
    paste("~/PhD/crows_project/results/BUSCO/run_",
          sp,
          "_busco/full_table_",
          sp,
          "_busco",sep=""),header=FALSE,fill = TRUE)
  colnames(dat)<-c("BUSCO","status","scaff_name","start","end","score","length")
  # Calculate the proportions in each status
  spres<-tapply(dat$BUSCO,INDEX = list(dat$status),length)/nrow(dat)
  res<-rbind(res,cbind(spres,sp,dimnames(spres)[[1]]))
  # Extract the "COMPLETE" ones from the table
  assign(paste(sp,"_compl",sep=""),dat[dat$status=="Complete",])
}
row.names(res)<-seq(1,nrow(res))
colnames(res)<-c("prop","sp","status")
res$prop<-as.numeric(as.character(res$prop))
res$sp<-factor(as.character(res$sp),
                  labels = c("C. cornix","C. corax","C. dauuricus",
                           "C. frugilegus","C. hawaiiensis","C. moneduloides",
                           "C. splendens","C. tasmanicus","C. brachyrhynchos",
                           "C. kubaryi","C. woodfordi","C. monedula","T. guttata"))

res$status<-as.character(res$status)
str(res)

# Plot the proportions
ggplot()+
  geom_bar(data=res,aes(sp,prop,fill=status),stat="identity")+
  xlab("")+
  ylab("Proportion")+
  scale_fill_manual("",values=c("#d7191c","#fdae61","#abd9e9","#2c7bb6"))+
  my.theme+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=10,face = "italic"),
        axis.text.y = element_text(size=10))

# Get common "Complete" genes
common<-Reduce(intersect,list(`3sp_Ccornix_compl`$BUSCO,
                      Ccorx_compl$BUSCO,Cdau_compl$BUSCO,
                      Cfru_compl$BUSCO, Chaw_compl$BUSCO,
                      Cmon_compl$BUSCO,Csple_compl$BUSCO,
                      Ctas_compl$BUSCO,Cbrach_compl$BUSCO,
                      Ckub_compl$BUSCO,Cwoo_compl$BUSCO,
                      Cmone_compl$BUSCO,Tgut_compl$BUSCO))
length(common)
# Write the common genes as a list
write.table(common,"~/PhD/crows_project/results/BUSCO/common_complete_genes.list",
            col.names=FALSE,row.names=FALSE,quote=FALSE)






