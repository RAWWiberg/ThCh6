###
# Plotting crow theta and FST results
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
options(scipen = 6)

or <-"#E69F00"
bl <-"#0072B2"
winsize<-"50kb"
anc<-""
# Read in data
####

## C. woodfordi
#####
Cwoo_dat<-read.table(paste("Cwoo_win",winsize,anc,"_thetas.tab",sep=""),
                     header=TRUE,sep="\t")
head(Cwoo_dat)
# Print as coordinate file
cwoodat_tab<-cbind(as.character(Cwoo_dat$Chr),
                   Cwoo_dat$WinCenter,
                   Cwoo_dat$WinCenter+1,
                   as.character(Cwoo_dat$Chr),
                   Cwoo_dat$WinCenter,
                   Cwoo_dat$nSites,
                   Cwoo_dat$tP,
                   Cwoo_dat$fayh,
                   Cwoo_dat$Tajima)
cwoodat_tab<-as.data.frame(cwoodat_tab)
head(cwoodat_tab)
write.table(cwoodat_tab,
            paste("~/PhD/crows_project/results/pop_genomics/Cwoo_win",
                  winsize,anc,".tab",sep=""),
            col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
head(Cwoo_dat)

# Read in the table file for windows with converted coordinates
Cwoo_dat<-read.table(paste("~/PhD/crows_project/results/pop_genomics/Cwoo_win",
                           winsize,anc,"_crowChr.tab",sep=""),
                     header=FALSE,sep="\t")
colnames(Cwoo_dat)<-c("chr","spos","epos","scaff","wincent",
                      "nsites","tP","fayH","tajD")
head(Cwoo_dat)
# Remove unlocalized and unplaced (Un) chromosomes
Cwoo_dat<-Cwoo_dat[grep("Un",Cwoo_dat$chr,invert=TRUE),]
# Sort by chromosome and by chromosome start
Cwoo_dat<-Cwoo_dat[order(Cwoo_dat$chr,Cwoo_dat$spos),]
Cwoo_dat$win<-paste(Cwoo_dat$scaff,"_",Cwoo_dat$wincent,sep="")
head(Cwoo_dat)
# Plot
length(unique(Cwoo_dat$chr))
unique(Cwoo_dat$chr)
# PLOT TAJIMA'S D
# Take only first 10 chromosomes
chrs<-head(unique(Cwoo_dat$chr),n=10)
chrs<-tail(unique(Cwoo_dat$chr),n=10)
# Order chromosomes
chrs<-c("1","1A","2","3","4","4A","5",
        "6","7","8","9","10",
        "11","12","13","14","15",
        "17","18","19","20",
        "21","22","23","24",
        "25","26","27","28","Z")
# Order scaffolds
scaffs<-as.character(unique(Cwoo_dat$scaff))
# Order chromosomes
Cwoo_dat$chr<-factor(Cwoo_dat$chr,levels=chrs)
Cwoo_dat$scaffold<-factor(Cwoo_dat$scaff,levels=scaffs)

str(Cwoo_dat)
head(Cwoo_dat[c(1,4,5)],n=200)
# PLOTTING: Cwoo
# ####
# Tajima's D
ggplot()+
  geom_line(data=Cwoo_dat,
            aes(epos/1000000,tajD,colour=scaff),
            size=0.1)+
  xlab("Position (Mb)")+
  ylab(expression(paste("Tajima's ",italic(D))))+
  scale_y_continuous(limits=c(-3,3),breaks=c(-3,-1,1,3))+  
  scale_colour_manual("",values=rep(c("black","red"),
                                    length.out=length(unique(Cwoo_dat$scaff))),
                      guide=FALSE)+
  geom_hline(yintercept=0,colour="grey")+
  facet_grid(chr~.)+
  my.theme+theme(axis.title = element_text(size=10),
                 axis.text.y = element_text(size=6),
                 axis.text.x = element_text(size=10),
                 strip.text.y = element_text(size=10,angle=0),
                 legend.text = element_text(size=8))
# Fay and Wu's H
ggplot()+
  geom_line(data=Cwoo_dat,
             aes(epos/1000000,fayH,colour=scaff),
             size=0.1)+
  xlab("Position (Mb)")+
  ylab(expression(paste("Fay and Wu's ",italic(H))))+
  scale_y_continuous(limits=c(-3.5,1.5),breaks=c(-3,-2,-1,0,1,2))+
  scale_colour_manual("",values=rep(c("black","red"),
                                    length.out=length(unique(Cwoo_dat$scaff))),
                      guide=FALSE)+
  geom_hline(yintercept=0,colour="grey")+
  facet_grid(chr~.)+
  my.theme+theme(axis.title = element_text(size=10),
                 axis.text.y = element_text(size=6),
                 axis.text.x = element_text(size=10),
                 strip.text.y = element_text(size=10,angle=0))
# Theta pi
ggplot()+
  geom_line(data=Cwoo_dat,
             aes(epos/1000000,tP/nsites,colour=scaff),
             size=0.1)+
  xlab("Position (Mb)")+
  ylab(expression(pi))+
  scale_y_continuous(limits=c(0,0.015),breaks=c(0,0.005,0.01))+
  scale_colour_manual("",values=rep(c("black","red"),
                                    length.out=length(unique(Cwoo_dat$scaff))),
                      guide=FALSE)+
  geom_hline(yintercept=0,colour="grey")+
  facet_grid(chr~.)+
  my.theme+theme(axis.title = element_text(size=10),
                 axis.text.y = element_text(size=6),
                 axis.text.x = element_text(size=10),
                 strip.text.y = element_text(size=10,angle=0))
# ####


## C. moneduloides
Cmon_dat<-read.table(paste("Cmon_win",winsize,anc,"_thetas.tab",sep=""),
                     header=TRUE,sep="\t")
# Print as coordinate file
cmondat_tab<-cbind(as.character(Cmon_dat$Chr),
                   Cmon_dat$WinCenter,
                   Cmon_dat$WinCenter+1,
                   as.character(Cmon_dat$Chr),
                   Cmon_dat$WinCenter,
                   Cmon_dat$nSites,
                   Cmon_dat$tP,
                   Cmon_dat$fayh,
                   Cmon_dat$Tajima)
cmondat_tab<-as.data.frame(cmondat_tab)
head(cmondat_tab)
write.table(cmondat_tab,
            paste("~/PhD/crows_project/results/pop_genomics/Cmon_win",
                  winsize,anc,".tab",sep=""),
            col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
head(Cmon_dat)

# Read in the table file for windows with converted coordinates
Cmon_dat<-read.table(paste("~/PhD/crows_project/results/pop_genomics/Cmon_win",
                           winsize,anc,"_crowChr.tab",sep=""),
                     header=FALSE,sep="\t")
colnames(Cmon_dat)<-c("chr","spos","epos","scaff","wincent",
                      "nsites","tP","fayH","tajD")
head(Cmon_dat)
# Remove unlocalized and unplaced (Un) chromosomes
Cmon_dat<-Cmon_dat[grep("Un",Cmon_dat$chr,invert=TRUE),]
# Sort by chromosome and by chromosome start
Cmon_dat<-Cmon_dat[order(Cmon_dat$chr,Cmon_dat$spos),]
Cmon_dat$win<-paste(Cmon_dat$scaff,"_",Cmon_dat$wincent,sep="")
head(Cmon_dat)
# Plot
length(unique(Cmon_dat$chr))
unique(Cmon_dat$chr)
# PLOT TAJIMA'S D
# Take only first 10 chromosomes
chrs<-head(unique(Cmon_dat$chr),n=10)
chrs<-tail(unique(Cmon_dat$chr),n=10)
# Order chromosomes
chrs<-c("1","1A","2","3","4","4A","5",
        "6","7","8","9","10",
        "11","12","13","14","15",
        "17","18","19","20",
        "21","22","23","24",
        "25","26","27","28","Z")
# Order scaffolds
scaffs<-as.character(unique(Cmon_dat$scaff))
# Order chromosomes
Cmon_dat$chr<-factor(Cmon_dat$chr,levels=chrs)
Cmon_dat$scaffold<-factor(Cmon_dat$scaff,levels=scaffs)

str(Cmon_dat)
head(Cmon_dat[c(1,4,5)],n=200)
# PLOTTING: Cmon
# ####
# Tajima's D
ggplot()+
  geom_line(data=Cmon_dat,
             aes(epos/1000000,tajD,colour=scaff),
             size=0.1)+
  xlab("Position (Mb)")+
  ylab(expression(paste("Tajima's ",italic(D))))+
  scale_y_continuous(limits=c(-3,3),breaks=c(-3,-1,1,3))+
  scale_colour_manual("",values=rep(c("black","red"),
                                    length.out=length(unique(Cmon_dat$scaff))),
                      guide=FALSE)+
  geom_hline(yintercept=0,colour="grey")+
  facet_grid(chr~.)+
  my.theme+theme(axis.title = element_text(size=10),
                 axis.text.y = element_text(size=6),
                 axis.text.x = element_text(size=10),
                 strip.text.y = element_text(size=10,angle=0))

# theta pi
ggplot()+
  geom_point(data=Cmon_dat,
            aes(epos/1000000,tP/nsites,colour=scaff),
            size=0.1)+
  xlab("Position (Mb)")+
  ylab(expression(pi))+
  scale_y_continuous(limits=c(0,0.015),breaks=c(0,0.005,0.01))+
  scale_colour_manual("",values=rep(c("black","red"),
                                    length.out=length(unique(Cmon_dat$scaff))),
                      guide=FALSE)+
  geom_hline(yintercept=0,colour="grey")+
  facet_grid(chr~.)+
  my.theme+theme(axis.title = element_text(size=10),
                 axis.text.y = element_text(size=6),
                 axis.text.x = element_text(size=10),
                 strip.text.y=element_text(size=10,angle=0))

# Fay and Wu's H
ggplot()+
  geom_point(data=Cmon_dat,
            aes(epos/1000000,fayH,colour=scaff),
            size=0.1)+
  xlab("Position (Mb)")+
  ylab(expression(paste("Fay and Wu's ",italic(H))))+
  scale_y_continuous(limits=c(-3.5,1.5),breaks=c(-3,-2,-1,0,1,2))+
  geom_hline(yintercept=0,colour="grey")+
  scale_colour_manual("",values=rep(c("black","red"),
                                    length.out=length(unique(Cmon_dat$scaff))),
                      guide=FALSE)+
  facet_grid(chr~.)+
  my.theme+theme(axis.title = element_text(size=10),
                axis.text.y = element_text(size=6),
                axis.text.x = element_text(size=10),
                strip.text.y=element_text(size=10))



####
# Both Species
#####

# Mean and median pi on each chromosome
# Combine Cwoo and Cmon
Cwoo_dat$sp<-rep("C. woodfordi",nrow(Cwoo_dat))
Cmon_dat$sp<-rep("C. moneduloides",nrow(Cmon_dat))
all_dat<-rbind(Cwoo_dat,Cmon_dat)

##
# PLOTTING BOTH
# ####
pi_m<-tapply(all_dat$tP/all_dat$nsites,INDEX=list(all_dat$chr,all_dat$sp),mean)
# Read in the chromosome sizes
crow_chr<-read.table("crowChr_sizes.tab",header=FALSE,sep="\t")
colnames(crow_chr)<-c("chr","length")
crow_chr<-crow_chr[crow_chr$chr %in% row.names(pi_m),]
#order the chromosomes
crow_chr<-crow_chr[order(crow_chr$chr),]
pi_m<-pi_m[order(row.names(pi_m)),]
nrow(crow_chr)
# Plot pi as a function of chromosome length
pi_dat<-data.frame(cwoopi=pi_m[,2],cmonpi=pi_m[,1],
                   length=crow_chr$length,chr=crow_chr$chr)
pi_dat_m<-melt(pi_dat,measure.vars = c("cwoopi","cmonpi"))
head(pi_dat_m)
ggplot()+
  geom_point(data=pi_dat_m[pi_dat_m$chr!="chr_Z",],
             aes(length/1000000,value,fill=variable),shape=21,size=3)+
  xlab("Length of Chromosome (Mb)")+
  ylab(expression(pi))+
  scale_fill_manual("",values=c(bl,or),
                      labels=c(expression(paste(bolditalic(C.),
                                                bolditalic( woodfordi))),
                               expression(paste(bolditalic(C.),
                                                bolditalic( moneduloides)))))+
  my.theme+theme(axis.title = element_text(size=11),
                 axis.text.y = element_text(size=11),
                 axis.text.x = element_text(size=11),
                 strip.text.y = element_text(size=11,angle=0),
                 legend.position = "top",
                 legend.text = element_text(size=11),
                 legend.key=element_rect(fill="white",colour="grey"))

cor.test(pi_dat$cwoopi[pi_dat$chr!="chr_Z" ],
         pi_dat$length[pi_dat$chr!="chr_Z"],method="spearman")
cor.test(pi_dat$cmonpi[pi_dat$chr!="chr_Z"],
         pi_dat$length[pi_dat$chr!="chr_Z"],method="spearman")
# ####

##
# WINDOW PLOTS: Both
# #####
# pi
ggplot()+
  geom_hline(yintercept=0,colour="grey")+
  geom_line(data=all_dat,
             aes(epos/1000000,tP/nsites,colour=sp),
             alpha=1/2)+
  xlab("Position (Mb)")+
  ylab(expression(pi))+
  scale_y_continuous(limits=c(0,0.015),breaks=c(0,0.005,0.01,0.015))+
  scale_colour_manual("",values=c("black","grey50"))+
  facet_grid(chr~.)+
  my.theme+theme(axis.title = element_text(size=10),
                 axis.text.y = element_text(size=6),
                 axis.text.x = element_text(size=10),
                 strip.text.y=element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.position = "top")
# Autosome
tapply(all_dat$tP[all_dat$chr!="Z"]/all_dat$nsites[all_dat$chr!="Z"],
       INDEX=list(all_dat$sp[all_dat$chr!="Z"]),mean)
tapply(all_dat$tP[all_dat$chr!="Z"]/all_dat$nsites[all_dat$chr!="Z"],
       INDEX=list(all_dat$sp[all_dat$chr!="Z"]),median)
tapply(all_dat$tP[all_dat$chr!="Z"]/all_dat$nsites[all_dat$chr!="Z"],
       INDEX=list(all_dat$sp[all_dat$chr!="Z"]),range)

# Z
tapply(all_dat$tP[all_dat$chr=="Z"]/all_dat$nsites[all_dat$chr=="Z"],
       INDEX=list(all_dat$sp[all_dat$chr=="Z"]),mean)
tapply(all_dat$tP[all_dat$chr=="Z"]/all_dat$nsites[all_dat$chr=="Z"],
       INDEX=list(all_dat$sp[all_dat$chr=="Z"]),median)
tapply(all_dat$tP[all_dat$chr=="Z"]/all_dat$nsites[all_dat$chr=="Z"],
       INDEX=list(all_dat$sp[all_dat$chr=="Z"]),range)

# Tajima's D
ggplot()+
  geom_hline(yintercept=0,colour="grey")+
  geom_line(data=all_dat,
            aes(epos/1000000,tajD,colour=sp),
            alpha=1/2)+
  xlab("Position (Mb)")+
  ylab(expression(paste("Tajima's ",italic(D))))+
  scale_colour_manual("",values=c("darkred","grey10"))+
  facet_grid(chr~.)+
  my.theme+theme(axis.title = element_text(size=10),
                 axis.text.y = element_text(size=6),
                 axis.text.x = element_text(size=10),
                 strip.text.y=element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.position = "top")
# Autosome
tapply(all_dat$tajD[all_dat$chr!="Z"],
       INDEX=list(all_dat$sp[all_dat$chr!="Z"]),mean)
tapply(all_dat$tajD[all_dat$chr!="Z"],
       INDEX=list(all_dat$sp[all_dat$chr!="Z"]),median)
tapply(all_dat$tajD[all_dat$chr!="Z"],
       INDEX=list(all_dat$sp[all_dat$chr!="Z"]),range)

# Z
tapply(all_dat$tajD[all_dat$chr=="Z"],
       INDEX=list(all_dat$sp[all_dat$chr=="Z"]),mean)
tapply(all_dat$tajD[all_dat$chr=="Z"],
       INDEX=list(all_dat$sp[all_dat$chr=="Z"]),median)
tapply(all_dat$tajD[all_dat$chr=="Z"],
       INDEX=list(all_dat$sp[all_dat$chr=="Z"]),range)

# Fay and Wu's H
ggplot()+
  geom_hline(yintercept=0,colour="grey")+
  geom_line(data=all_dat,
            aes(epos/1000000,fayH,colour=sp),
            alpha=1/2)+
  xlab("Position (Mb)")+
  ylab(expression(paste("Fay and Wu's ",italic(H))))+
  scale_colour_manual("",values=c("red","grey10"))+
  facet_grid(chr~.)+
  my.theme+theme(axis.title = element_text(size=10),
                 axis.text.y = element_text(size=6),
                 axis.text.x = element_text(size=10),
                 strip.text.y=element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.position = "top")
# Autosome
tapply(all_dat$fayH[all_dat$chr!="Z"],
       INDEX=list(all_dat$sp[all_dat$chr!="Z"]),mean)
tapply(all_dat$fayH[all_dat$chr!="Z"],
       INDEX=list(all_dat$sp[all_dat$chr!="Z"]),median)
tapply(all_dat$fayH[all_dat$chr!="Z"],
       INDEX=list(all_dat$sp[all_dat$chr!="Z"]),range)

# Z
tapply(all_dat$fayH[all_dat$chr=="Z"],
       INDEX=list(all_dat$sp[all_dat$chr=="Z"]),mean)
tapply(all_dat$fayH[all_dat$chr=="Z"],
       INDEX=list(all_dat$sp[all_dat$chr=="Z"]),median)
tapply(all_dat$fayH[all_dat$chr=="Z"],
       INDEX=list(all_dat$sp[all_dat$chr=="Z"]),range)

# ####


# FST: Cmon Cwoo
# #####
# Read in the FST data
CmonCwoo_dat<-read.table("Cmon_Cwoo_Fst_slidingwindow_red.tab",
                     header=TRUE,sep="\t")
head(CmonCwoo_dat)
# Remove chrM
CmonCwoo_dat<-CmonCwoo_dat[grep("chrM",CmonCwoo_dat$chr,invert=TRUE),]
# Print as coordinate file
cmoncwoofst_tab<-cbind(as.character(CmonCwoo_dat$chr),
                       CmonCwoo_dat$midPos,
                       CmonCwoo_dat$midPos+1,
                   as.character(CmonCwoo_dat$chr),
                   CmonCwoo_dat$midPos,
                   CmonCwoo_dat$Nsites,
                   CmonCwoo_dat$FST)
cmoncwoofst_tab<-as.data.frame(cmoncwoofst_tab)
head(cmoncwoofst_tab)
write.table(cmoncwoofst_tab,
            "Cmon_Cwoo_Fst_sliwin.tab",
            col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

# Read in the gtf file for windows with converted coordinates
CmonCwoo_dat<-read.table(
  "~/PhD/crows_project/results/pop_genomics/Cmon_Cwoo_Fst_win50kb_crowChr.tab",
                     header=FALSE,sep="\t")
colnames(CmonCwoo_dat)<-c("chr","spos","epos","scaff","wincent","nsites","fst")
head(CmonCwoo_dat)
# Remove unlocalize, unplaced (Up)
CmonCwoo_dat<-CmonCwoo_dat[grep("Un",CmonCwoo_dat$chr,invert=TRUE),]
CmonCwoo_dat$win<-paste(CmonCwoo_dat$scaff,"_",CmonCwoo_dat$wincent,sep="")

# Sort by chromosome and by chromosome start
CmonCwoo_dat<-CmonCwoo_dat[order(CmonCwoo_dat$chr,CmonCwoo_dat$spos),]
# remove the "linkage group" (LG) scaffold
unique(CmonCwoo_dat$chr)
# Order chromosomes
chrs<-c("1","1A","2","3","4","4A","5",
        "6","7","8","9","10",
        "11","12","13","14","15",
        "17","18","19","20",
        "21","22","23","24",
        "25","26","27","28","Z")

# Order chromosomes and scaffolds
CmonCwoo_dat$chr<-factor(CmonCwoo_dat$chr,levels=chrs)
CmonCwoo_dat$scaff<-factor(CmonCwoo_dat$scaff,
                           levels=as.character(unique(CmonCwoo_dat$scaff)))
# ####

# PLOTTING: FST
# ####
ggplot()+
  geom_line(data=CmonCwoo_dat,
             aes((epos-100)/1000000,fst,colour=scaff))+
  xlab("Position (Mb)")+
  ylab(expression(F[st]))+
  scale_colour_manual("",values=rep(c("black","grey50"),
                                    length.out=length(unique(CmonCwoo_dat$scaff))),
                      guide=FALSE)+
  facet_grid(chr~.)+
  my.theme+theme(axis.title = element_text(size=10),
                 axis.text.y = element_text(size=6),
                 axis.text.x = element_text(size=10),
                 strip.text.y = element_text(size=10,angle=0))
ggplot()+
  geom_boxplot(data=CmonCwoo_dat,
               aes(chr,fst),outlier.size = 0.3)+
  xlab(expression(bold("Chromosome")))+
  ylab(expression(bolditalic(F)[bold(st)]))+
  my.theme+theme(axis.title = element_text(size=11),
                 axis.text.y = element_text(size=11),
                 axis.text.x = element_text(size=11,angle=45,hjust=1,vjust=1),
                 strip.text.y = element_text(size=11,angle=0),
                 legend.text = element_text(size=11))
  
tapply(CmonCwoo_dat$fst,INDEX = list(CmonCwoo_dat$chr),mean)
fst_m<-tapply(CmonCwoo_dat$fst,INDEX = list(CmonCwoo_dat$chr),median)
summary(fst_m[seq(1,nrow(fst_m)-1)])
fst_m[nrow(fst_m)]
nrow(CmonCwoo_dat)
head(CmonCwoo_dat)
# ####

# PLOT FST AGAINST PI
# ####
head(CmonCwoo_dat)
nrow(CmonCwoo_dat)
head(Cmon_dat)
nrow(Cmon_dat)
# Subset
CmonCwoo_dat_fstpi<-CmonCwoo_dat[which((CmonCwoo_dat$win %in% Cmon_dat$win)),]
#Cmon
Cmon_dat_fstpi<-Cmon_dat[which((Cmon_dat$win %in% CmonCwoo_dat$win)),]
nrow(Cmon_dat_fstpi)
head(Cmon_dat_fstpi)
#Cwoo
Cwoo_dat_fstpi<-Cwoo_dat[which((Cwoo_dat$win %in% CmonCwoo_dat$win)),]
nrow(Cwoo_dat_fstpi)
head(Cwoo_dat_fstpi)
#Combine
Cmon_dat_fstpi$fst<-CmonCwoo_dat_fstpi$fst
Cwoo_dat_fstpi$fst<-CmonCwoo_dat_fstpi$fst
Cwoo_dat_fstpi$sp<-rep("C. woodfordi",nrow(Cwoo_dat_fstpi))
Cmon_dat_fstpi$sp<-rep("C. moneduloides",nrow(Cmon_dat_fstpi))
all_dat_fstpi<-rbind(Cwoo_dat_fstpi,Cmon_dat_fstpi)
str(all_dat_fstpi)
all_dat_fstpi$sp<-factor(all_dat_fstpi$sp)
unique(all_dat_fstpi$sp)
head(all_dat_fstpi)
all_dat_fstpi$sp<-factor(all_dat_fstpi$sp,
                         labels=c(expression(paste(bolditalic(C.),
                                                   bolditalic( moneduloides))),
                                  expression(paste(bolditalic(C.),
                                                   bolditalic( woodfordi)))))
str(all_dat_fstpi)
unique(all_dat_fstpi$sp)
# PLot
ggplot()+geom_point(data=all_dat_fstpi[all_dat_fstpi$chr!="Z",],
  aes(x=fst,y=tP/nsites,fill=sp),colour="black",
  alpha=1/4,shape=21)+
  scale_fill_manual("",labels=c(expression(paste(italic(C.),italic( moneduloides))),
                                expression(paste(italic(C.),italic( woodfordii)))),
                    values=c(or,bl),guide=FALSE)+
  xlab(expression(bolditalic(F)[bold(st)]))+
  ylab(expression(bold("\u03C0")))+
  facet_grid(sp~.,labeller = label_parsed)+
  my.theme+theme(axis.title = element_text(size=11,face="bold"),
                 axis.text.y = element_text(size=11),
                 axis.text.x = element_text(size=11),
                 strip.text.y = element_text(size=11,angle=270,face="bold"),
                 legend.text = element_text(size=11),
                 legend.position="top")

# Plot pi for orhologous windows in cmon and cwoo
plot((Cmon_dat_fstpi$tP/Cmon_dat_fstpi$nsites),
     (Cwoo_dat_fstpi$tP/Cwoo_dat_fstpi$nsites))

# ####

# Combine TajD, pi and FST datasets for a multipanel plot

head(all_dat)
str(all_dat)

all_dat$pi<-all_dat$tP/all_dat$nsites
# Melt all_dat
all_dat_m<-melt(data = all_dat,measure.vars = c("pi","fayH","tajD"))
head(all_dat_m)
all_dat_m$variable<-factor(all_dat_m$variable,
                           labels=c(expression(pi),
                           expression(paste("Fay & Wu's ",italic(H))),
                           expression(paste("Tajima's ",italic(D)))))
ggplot()+
  geom_boxplot(data=all_dat_m,aes(chr,value,fill=sp),
               outlier.size = 0.2)+
  facet_grid(variable~.,scales = "free_y",labeller = label_parsed)+
  ylab("")+
  xlab("Chromosome")+
  scale_fill_manual("",
                    labels=c(expression(paste(bolditalic(C.),
                                              bolditalic( moneduloides))),
                               expression(paste(bolditalic(C.),
                                                bolditalic( woodfordi)))),
                    values=c(or,bl))+
  my.theme+theme(axis.title = element_text(size=10,face="bold"),
                 axis.text.y = element_text(size=10),
                 axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
                 strip.text.y = element_text(size=10,angle=270,face="bold"),
                 legend.text = element_text(size=10),
                 legend.position="top",
                 legend.key=element_rect(fill="white",colour="grey"))




###
# DO SOME SIMULATIONS
# #####
# Read in the data
nsnps<-50000
nsims<-100
sim_dat<-read.table(
  paste("~/PhD/crows_project/results/pop_genomics/crows_pop_gen_sims_ns",nsnps,"_s",
        nsims,".tab",sep=""),header=TRUE,sep="\t")

# Re-code the factors for nice labels
sim_dat$N<-factor(sim_dat$N,
                  levels=paste("N = ",c(5,10,50,100),sep=""),
                  labels=c(expression(paste(italic(N)," = 5")),
                           expression(paste(italic(N)," = 10")),
                           expression(paste(italic(N)," = 50")),
                           expression(paste(italic(N)," = 100"))))
sim_dat$Ne<-factor(sim_dat$Ne,levels=c("Ne = 1000",
                                       "Ne = 10000",
                                       "Ne = 20000",
                                       "Ne = 100000",
                                       "Ne = 200000",
                                       "Ne = 1000000",
                                       "Ne = 2000000"),
                   labels = c(expression(paste(N[e],"= 1,000")),
                              expression(paste(N[e],"= 10,000")),
                              expression(paste(N[e],"= 20,000")),
                              expression(paste(N[e],"= 100,000")),
                              expression(paste(N[e],"= 200,000")),
                              expression(paste(N[e],"= 1,000,000")),
                              expression(paste(N[e],"= 2,000,000"))))
# Make mu and alpha vectors into factors
sim_dat$muf<-factor(sim_dat$mu,levels=c("0.00000000123",
                                        "0.00000000191",
                                        "0.00000000221",
                                        "0.0000000023"),
                    labels = c(expression(paste("1.23x",10^-9)),
                               expression(paste("1.91x",10^-9)),
                               expression(paste("2.21x",10^-9)),
                               expression(paste("2.3x",10^-9))))
sim_dat$alphaf<-factor(sim_dat$alpha)

# PLOT TAJIMA'S D FOR EACH PARAMETER COMBINATION
head(sim_dat)
ggplot()+
  geom_boxplot(data=sim_dat,aes(N,tajd,fill=muf))+
  facet_grid(Ne~.,labeller = label_parsed)+
  ylab(expression(paste("Tajima's ",italic(D))))+
  ylim(-2,2)+
  scale_fill_manual("Mutation Rate",
                    values=c("darkred","red","yellow","white"),
                      labels=c(expression(paste("1.23x",10^-9)),
                               expression(paste("1.91x",10^-9)),
                               expression(paste("2.21x",10^-9)),
                               expression(paste("2.3x",10^-9))))+
  xlab("")+
  scale_x_discrete(labels=c(expression(paste(italic(N)," = 5")),
                            expression(paste(italic(N)," = 10")),
                            expression(paste(italic(N)," = 50")),
                            expression(paste(italic(N)," = 100"))))+
  my.theme+theme(axis.title = element_text(size=10,face="bold"),
                 axis.text.y = element_text(size=10),
                 axis.text.x = element_text(size=10),
                 strip.text.y=element_text(angle=0,size=10),
                 strip.text.x=element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.title = element_text(size=10,face="bold"),
                 legend.position = "top")


# PLOT ALPHA VS. PI FOR EACH PARAMETER COMBINATION
ggplot()+
  geom_point(data=sim_dat,aes(alpha,pi,fill=muf),shape=21)+
  facet_grid(Ne~N,labeller = label_parsed)+
  xlab(expression(paste("4",N[e],mu,sep="")))+
  ylab(expression(pi))+
  scale_y_continuous(breaks=c(0, 0.010, 0.020))+
  scale_fill_manual("Mutation Rate",
                    values=c("darkred","red","yellow","white"),
                    labels=c(expression(paste("1.23x",10^-9)),
                             expression(paste("1.91x",10^-9)),
                             expression(paste("2.21x",10^-9)),
                             expression(paste("2.3x",10^-9))))+
  my.theme+theme(axis.title = element_text(size=10,face="bold"),
                 axis.text.y = element_text(size=10),
                 axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=10),
                 strip.text.y=element_text(angle=0,size=10),
                 strip.text.x=element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.title = element_text(size=10,face="bold"),
                 legend.position = "top")


tapply(sim_dat$tajd,INDEX = list(sim_dat$Ne,sim_dat$N,sim_dat$mu),mean,na.rm=TRUE)
tapply(sim_dat$pi,INDEX = list(sim_dat$Ne,sim_dat$N,sim_dat$mu),mean,na.rm=TRUE)

# Get the 1st and 5th percentiles of the Tajima's D distribution

q1tajd<-tapply(sim_dat$tajd,INDEX = list(sim_dat$mu,sim_dat$N,sim_dat$Ne),quantile,
       probs=c(0.01), na.rm=TRUE)
dim(q1tajd)
dimnames(q1tajd)
# Get the range of the 1st percentile
for(i in 1:dim(q1tajd)[3]){
  ne<-dimnames(q1tajd)[[3]][i]
  for(j in 1:dim(q1tajd)[2]){
    n<-dimnames(q1tajd)[[2]][j]
    q1<-q1tajd[,j,i]
    cat(ne,n,range(q1),"\n")
  }
}

tapply(sim_dat$tajd,INDEX = list(sim_dat$Ne,sim_dat$N,sim_dat$mu),quantile,
       probs=c(0.05), na.rm=TRUE)
# ####





###
# Get regions of low Tajima's D
###
# How many windows show Tajima's D values < tajd_thresh
tajd_thresh<--1.6
nrow(Cmon_dat[Cmon_dat$tajD < tajd_thresh,])

head(Cmon_dat)
head(Cwoo_dat)
head(CmonCwoo_dat)

#Extract the region +- bMb
wins<-50000
b<-1000000
TajD_wins<-Cmon_dat[Cmon_dat$tajD < tajd_thresh,]
# Expand the start and end positions of windows to match their width
TajD_wins$spos<-TajD_wins$spos-1000000
TajD_wins$epos<-TajD_wins$epos+1000000
head(TajD_wins)

# Write the location as a bed format table.
write.table(TajD_wins[,c(1,2,3)],
            paste(
              "~/PhD/crows_project/results/pop_genomics/Cmon_lowTajD_wins_",
              winsize,".bed",sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


# Load the merged bed table.
TajD_wins<-read.table(
  "~/PhD/crows_project/results/pop_genomics/Cmon_lowTajD_wins_merged.bed",
  sep="\t")
colnames(TajD_wins)<-c("chr","spos","epos")
TajD_wins$win<-paste(TajD_wins$chr,":",TajD_wins$spos,"-",TajD_wins$epos,sep="")
TajD_wins
head(Cmon_dat)
# For each region, get the windows that fall within the interval.
TajD_regions<-data.frame()
CmonTajD_regions<-data.frame()
CwooTajD_regions<-data.frame()
CmonCwoo_fst_regions<-data.frame()
for(win in TajD_wins$win){
  chr<-as.character(TajD_wins$chr[TajD_wins$win==win])
  spos<-TajD_wins$spos[TajD_wins$win==win]
  epos<-TajD_wins$epos[TajD_wins$win==win]
  # Get regions from Cmon
  cmonregions<-Cmon_dat[Cmon_dat$chr==chr &
                         Cmon_dat$spos > spos &
                         Cmon_dat$epos < epos,]
  #Add an abspos and win id column
  cmonregions$abspos<-c(0,cumsum(diff(cmonregions$spos)))
  cmonregions$win<-rep(win,nrow(cmonregions))
  # Get regions from Cwoo
  cwooregions<-Cwoo_dat[Cwoo_dat$chr==chr &
                         Cwoo_dat$spos > spos &
                         Cwoo_dat$spos < epos,]
  #Add an abspos and win id column
  cwooregions$abspos<-c(0,cumsum(diff(cwooregions$spos)))
  cwooregions$win<-rep(win,nrow(cwooregions))
  # Get FST regions
  fst_regions<-CmonCwoo_dat[CmonCwoo_dat$chr==chr & 
                              CmonCwoo_dat$spos > spos &
                              CmonCwoo_dat$epos < epos,]
  fst_regions$abspos<-c(0,cumsum(diff(fst_regions$spos)))
  fst_regions$win<-rep(win,nrow(fst_regions))
  # Combine the data
  TajD_regions<-rbind(TajD_regions,cmonregions,cwooregions)
  CmonCwoo_fst_regions<-rbind(CmonCwoo_fst_regions,fst_regions)
}
#Organise the regions in order
TajD_regions$win<-factor(TajD_regions$win,
                         levels=c("1A:52349998-54399999","2:82255643-84315644",
                                  "2:138372818-140422819","4:30139573-32189574",
                                  "4:50709056-52759057","5:24519944-26569945",
                                  "5:26829944-29089945","5:46499944-48549945",
                                  "23:1518721-3578722","Z:4474020-6554021",
                                  "Z:71555854-73605855"))
TajD_regions
unique(TajD_regions$win)

# Plot: Tajima's D
ggplot()+
  geom_hline(yintercept = 0,colour="grey")+
  geom_line(data=TajD_regions,
            aes(abspos,tajD,colour=sp),size=0.8)+
  scale_colour_manual("",values=c(or,bl))+
  ylim(-2,2)+
  facet_grid(win~.)+
  xlab("")+
  ylab(expression(paste("Tajima's ",italic(D))))+
  my.theme+theme(axis.title = element_text(size=10,face="bold"),
                 axis.text.y = element_text(size=10),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 strip.text.y=element_text(angle=0,size=10),
                 strip.text.x=element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.title = element_text(size=10,face="bold"),
                 legend.position = "top")

# Plot: Fay and Wu's H
ggplot()+
  geom_hline(yintercept = 0,colour="grey")+
  geom_line(data=TajD_regions,
            aes(abspos,fayH,colour=sp))+
  scale_colour_manual("",values=c(or,bl))+
  facet_grid(win~.)+
  xlab("")+
  ylab(expression(paste("Fay and Wu's ",italic(H))))+
  my.theme+theme(axis.title = element_text(size=10,face="bold"),
                 axis.text.y = element_text(size=10),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 strip.text.y=element_text(angle=0,size=10),
                 strip.text.x=element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.title = element_text(size=10,face="bold"),
                 legend.position = "top")
# Plot FST
head(CmonCwoo_fst_regions)
ggplot()+
  geom_hline(yintercept = 0,colour="grey")+
  geom_line(data=CmonCwoo_fst_regions,
            aes(abspos,fst))+
  facet_grid(win~.)+
  xlab("")+
  ylab(expression(F[st]))+
  my.theme+theme(axis.title = element_text(size=10,face="bold"),
                 axis.text.y = element_text(size=10),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 strip.text.y=element_text(angle=0,size=10),
                 strip.text.x=element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.title = element_text(size=10,face="bold"),
                 legend.position = "top")

# Load the overlapping genes.
# ALL GENES
genes<-read.table(
  "~/PhD/crows_project/results/pop_genomics/Cmon_lowTajD_wins_merged_genesOverlap.bed",
                  header=FALSE,sep="\t")
colnames(genes)<-c("chr","start","end","id","V","strand","source","type","V2","name",
                   "winchr","winspos","winepos","nrbases")
genes$win<-paste(genes$winchr,":",genes$winspos,"-",genes$winepos,sep="")
head(genes)
# ZEBRA GENES
zebr_genes<-read.table(
  "~/PhD/crows_project/results/pop_genomics/Cmon_lowTajD_wins_merged_genesOverlap_ZEBRA.tab",
  header=TRUE,sep="\t")
head(zebr_genes)
zebr_genes<-zebr_genes[
  as.character(zebr_genes$Name) == as.character(zebr_genes$ensembl_ID),]

#Beak genes
beak_genes<-c("DKK2","CALM1","FOXO6","CHST11","TRPS1","NPR2")

head(TajD_wins)
head(TajD_regions)
# Give new positions to the genes
genes$regstart<-vector(length=nrow(genes))
genes$regend<-vector(length=nrow(genes))
for(win in TajD_wins$win){
  for(gene in genes$name[genes$win==win]){
    # Get the first position of the region
    start<-min(TajD_regions$spos[TajD_regions$win==win])
    # Calculate the new start position for each gene
    new_start<-genes$start[genes$win==win]-start
    g_length<-genes$end[genes$win==win]-genes$start[genes$win==win]
    new_end<-new_start+g_length

    # If any genes have a new start < 0
    new_start[which(new_start < 0)] <- 0
    # If any genes have a new end > max
    new_end[which(new_start > max(
      TajD_regions$abspos[TajD_regions$win==win]))]<-max(
        TajD_regions$abspos[TajD_regions$win==win])
    genes$regstart[genes$win==win]<-new_start
    genes$regend[genes$win==win]<-new_end
  }
}
head(genes)
# Reorder the factor levels
genes$win<-factor(genes$win,
                  levels=c("1A:52349998-54399999","2:82255643-84315644",
                           "2:138372818-140422819","4:30139573-32189574",
                           "4:50709056-52759057","5:24519944-26569945",
                           "5:26829944-29089945","5:46499944-48549945",
                           "23:1518721-3578722","Z:4474020-6554021",
                           "Z:71555854-73605855"))
# Remove anything with "LOC
nrow(genes[grep("LOC",genes$name,invert=TRUE),])
tgenes<-genes[grep("LOC",genes$name,invert=TRUE),]
# Subset ALL genes to thos in Zebra
genes[which(genes$name %in% zebr_genes$Name),]
head(TajD_regions)
# Tajima's D
ggplot()+
  geom_hline(yintercept = 0,colour="grey")+
  geom_line(data=TajD_regions,
            aes(abspos,tajD,colour=sp),size=0.8)+
  scale_colour_manual("",values=c(or,bl),
                      labels=c(expression(paste(bolditalic("C. "),
                                                bolditalic("moneduloides"))),
                               expression(paste(bolditalic("C. "),
                                                bolditalic("woodfordi")))))+
  geom_rect(data=genes[which(genes$name %in% zebr_genes$Name),],
            aes(ymin=1.8,ymax=2,xmin=regstart,xmax=regend))+
  geom_text(data=genes[which(genes$name %in% zebr_genes$Name),],
            aes(y=2.2,x=regstart,label=name),
            angle=0,size=2.2)+
  facet_grid(win~.)+
  xlab("")+
  ylim(-2,2.4)+
  ylab(expression(paste(bold("Tajima's "),bolditalic(D))))+
  my.theme+theme(axis.title = element_text(size=10,face="bold"),
                 axis.text.y = element_text(size=10),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 strip.text.y=element_text(angle=0,size=10),
                 strip.text.x=element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.title = element_text(size=10,face="bold"),
                 legend.position = "top",
                 legend.key = element_rect(fill="white",colour="grey"),
                 plot.title=element_text(size=10))

# Fay and Wu's H
ggplot()+
  geom_hline(yintercept = 0,colour="grey")+
  geom_line(data=TajD_regions,
            aes(abspos,fayH,colour=sp),size=0.8)+
  scale_colour_manual("",values=c(or,bl),
                      labels=c(expression(paste(bolditalic("C. "),
                                                bolditalic("moneduloides"))),
                               expression(paste(bolditalic("C. "),
                                                bolditalic("woodfordi")))))+
  geom_rect(data=genes[which(genes$name %in% beak_genes),],
            aes(ymin=0.8,ymax=1,xmin=regstart,xmax=regend))+
  geom_text(data=genes[which(genes$name %in% beak_genes),],
            aes(y=1.4,x=regstart,label=name),
            angle=0,size=2.2)+
  facet_grid(win~.)+
  ylim(-4,2)+
  xlab("")+
  ylab(expression(paste(bold("Fay and Wu's "),bolditalic(H))))+
  my.theme+theme(axis.title = element_text(size=10,face="bold"),
                 axis.text.y = element_text(size=10),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 strip.text.y=element_text(angle=0,size=10),
                 strip.text.x=element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.title = element_text(size=10,face="bold"),
                 legend.margin = margin(t=0.05,r=0.05,b=0.05,l=0.05),
                 legend.position = "top",
                 legend.key = element_rect(fill="white",colour="grey"),
                 plot.title=element_text(size=10))


