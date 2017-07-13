###
# Assessing the crow ortholog sets
# from orthAgogue output
#
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

# Load data: All species #####
gene_dat <- read.table("~/RData/all_Corvus_summary_stats.fa_summary",header = TRUE,
                       sep = "\t")
head(gene_dat)
# ####

# Analysis: All species ####
# %completeness = (1-Ns/length)*100
gene_dat$compl <- (1-(gene_dat$N_count/gene_dat$length))*100
head(gene_dat)
#total genes for each species
tapply(gene_dat$compl,gene_dat$species, length)

compl_plot <- ggplot()+
  geom_histogram(data = gene_dat,aes(compl,group= species))+
  facet_grid(species~.)
compl_plot+my.theme+theme(strip.text.y = element_text(angle = 0))

str(gene_dat)
#total counts for each species
aggregate(gene_dat[,c(3,4,5,6,7,8,9,10)],by=list(gene_dat$species),sum)

#number of genes with > 80% completeness
compl20<-melt(tapply(gene_dat$compl[gene_dat$compl > 20], 
       gene_dat$species[gene_dat$compl > 20], length))
colnames(compl20) <- c("species","count")
compl20$compl<-rep("20%",nrow(compl20))
#sequencing effort = N .bam files
compl20$seqeff <- c(15,32,1,16,16,1,1,11,10,10,1)

compl50<-melt(tapply(gene_dat$compl[gene_dat$compl > 50], 
                     gene_dat$species[gene_dat$compl > 50], length))
colnames(compl50) <- c("species","count")
compl50$compl<-rep("50%",nrow(compl50))
#sequencing effore = N .bam files
compl50$seqeff <- c(15,32,1,16,16,1,1,11,10,10,1)

compl90<-melt(tapply(gene_dat$compl[gene_dat$compl > 90], 
                     gene_dat$species[gene_dat$compl > 90], length))
colnames(compl90) <- c("species","count")
compl90$compl<-rep("90%",nrow(compl95))
#sequencing effore = N .bam files
compl95$seqeff <- c(15,32,1,16,16,1,1,11,10,10,1)

compldat <- rbind(compl20,compl50,compl90)
head(compldat)
#completeness plot
compl_plot <- ggplot() +
  geom_bar(data = compldat, aes(species,count,
                                group = compl),
           stat="identity",width = 0.5)+
  geom_text(data=compldat,aes(x=species,y=count+2000,
                              label=as.character(count)),size = 3)+
#  scale_y_continuous(limits = c(0,30000))+
  xlab("")+
  ylab("Count")+
  facet_grid(compl~.,scales = "free_y")
compl_plot+my.theme+theme(strip.text.y = element_text(angle = 0),
                          axis.text.x = element_text(size = 10,angle = 45,
                                                     vjust = 1,hjust=1),
                          axis.text.y = element_text(size = 14))
# ####

# How many genes for PAML analyses ####
# Tree1: 5 species
# Ccorone3sp, Cdau, Cmon, Cfru, Csple, Tgut
specs5 <- c("3sp-Ccornix","Cdau","Cfru","Csple","Cmon")
sp5_compldat<-gene_dat[gene_dat$species %in% specs5 & gene_dat$compl > 90,
                       c(1,2)]
nrow(sp5_compldat)
head(sp5_compldat)
#extract individual lists
sp3<-sp5_compldat$sequence_id[sp5_compldat$species == "3sp-Ccornix"]
Cdau<-sp5_compldat$sequence_id[sp5_compldat$species == "Cdau"]
Cfru<-sp5_compldat$sequence_id[sp5_compldat$species == "Cfru"]
Csple<-sp5_compldat$sequence_id[sp5_compldat$species == "Csple"]
Cmon<-sp5_compldat$sequence_id[sp5_compldat$species == "Cmon"]
#check overlap between lists
common<-Reduce(intersect,list(sp3,Cdau,Cfru,Csple,Cmon))
length(common)

#print "common" to a file
#write.table(common,
#          file = "~/PhD/crows_project/data/genes/5sp_common_genes.txt",
#          quote = FALSE,
#          col.names = FALSE,
#          row.names = FALSE)

# Tree2: 6 species
# Ccorone3sp, Cdau, Cmon, Cfru, Csple,Ctas, Tgut 
specs6 <- c("3sp-Ccornix","Cdau","Cfru","Ctas","Csple","Cmon")
sp6_compldat<-gene_dat[gene_dat$species %in% specs6& gene_dat$compl > 90,
                       c(1,2)]
nrow(sp6_compldat)
#extract individual lists
sp3<-sp6_compldat$sequence_id[sp6_compldat$species == "3sp-Ccornix"]
Cdau<-sp6_compldat$sequence_id[sp6_compldat$species == "Cdau"]
Cfru<-sp6_compldat$sequence_id[sp6_compldat$species == "Cfru"]
Csple<-sp6_compldat$sequence_id[sp6_compldat$species == "Csple"]
Cmon<-sp6_compldat$sequence_id[sp6_compldat$species == "Cmon"]
Ctas<-sp6_compldat$sequence_id[sp6_compldat$species == "Ctas"]
#check overlap between lists
common<-Reduce(intersect,list(sp3,Cdau,Cfru,Csple,Cmon,Ctas))
length(common)

#print "common" to a file
#write.table(common,
#            file = "~/PhD/crows_project/data/genes/6sp_common_genes.txt",
#            quote = FALSE,
#            col.names = FALSE,
#            row.names = FALSE)


# Tree3: 7 species
# Ccorone3sp, Cdau, Cmon, Cfru, Csple, Ctas, Ccorx, Tgut 
specs7 <- c("3sp-Ccornix","Cdau","Cfru","Ctas","Ccorx","Csple","Cmon")
sp7_compldat<-gene_dat[gene_dat$species %in% specs7 & gene_dat$compl > 90,
                       c(1,2)]
nrow(sp7_compldat)

#extract individual lists
sp3<-sp7_compldat$sequence_id[sp7_compldat$species == "3sp-Ccornix"]
Cdau<-sp7_compldat$sequence_id[sp7_compldat$species == "Cdau"]
Cfru<-sp7_compldat$sequence_id[sp7_compldat$species == "Cfru"]
Csple<-sp7_compldat$sequence_id[sp7_compldat$species == "Csple"]
Cmon<-sp7_compldat$sequence_id[sp7_compldat$species == "Cmon"]
Ctas<-sp7_compldat$sequence_id[sp7_compldat$species == "Ctas"]
Ccorx<-sp7_compldat$sequence_id[sp7_compldat$species == "Ccorx"]
#check overlap between lists
common<-Reduce(intersect,list(sp3,Cdau,Cfru,Csple,Cmon,Ctas,Ccorx))
length(common)

#print "common" to a file
#write.table(common,
#            file = "~/PhD/crows_project/data/genes/7sp_common_genes.txt",
#            quote = FALSE,
#            col.names = FALSE,
#            row.names = FALSE)

# Tree3: 8 species
# Ccorone3sp, Cdau, Cmon, Cfru, Csple, Ctas, Ccorx, Tgut 
specs8 <- c("3sp-Ccornix","Chaw","Cdau","Cfru","Ctas","Ccorx","Csple","Cmon")
sp8_compldat<-gene_dat[gene_dat$species %in% specs8 & gene_dat$compl > 90,
                       c(1,2)]
nrow(sp8_compldat)

#extract individual lists
sp3<-sp8_compldat$sequence_id[sp8_compldat$species == "3sp-Ccornix"]
Cdau<-sp8_compldat$sequence_id[sp8_compldat$species == "Cdau"]
Cfru<-sp8_compldat$sequence_id[sp8_compldat$species == "Cfru"]
Csple<-sp8_compldat$sequence_id[sp8_compldat$species == "Csple"]
Cmon<-sp8_compldat$sequence_id[sp8_compldat$species == "Cmon"]
Chaw<-sp8_compldat$sequence_id[sp8_compldat$species == "Chaw"]
Ctas<-sp8_compldat$sequence_id[sp8_compldat$species == "Ctas"]
Ccorx<-sp8_compldat$sequence_id[sp8_compldat$species == "Ccorx"]
#check overlap between lists
common<-Reduce(intersect,list(sp3,Cdau,Cfru,Csple,Cmon,Ctas,Ccorx,Chaw))
length(common)

#print "common" to a file
#write.table(common,
#            file = "~/PhD/crows_project/data/genes/8sp_common_genes.txt",
#            quote = FALSE,
#            col.names = FALSE,
#            row.names = FALSE)
# ####

#
# Ccornix v. Tguttata: For PAML Analysis
#
# Checking orthAgogue ortholog assignment
# ####
orth_dat <- read.table("~/RData/ccornix_v_tgut_orthAg_orthologs2.tab",
                       header = FALSE,sep = "\t")
colnames(orth_dat)<-c("sp1","sp2","score")
head(orth_dat)

# Nr Ccornix genes
length(unique(orth_dat$sp1[grep("Ccornix",orth_dat$sp1)]))

length(unique(orth_dat$sp2[grep("Ccornix",orth_dat$sp2)]))

head(orth_dat$sp1[grep("Ccornix",orth_dat$sp1)])

# Subset rows
ccornix_sp1<-orth_dat$sp1[grep("Ccornix",orth_dat$sp1)]
ccornix_sp2<-orth_dat$sp2[grep("Ccornix",orth_dat$sp2)]
length(unique(ccornix_sp1))
length(unique(ccornix_sp2))


# Nr Tgut genes
length(unique(orth_dat$sp1[grep("Tgut",orth_dat$sp1)]))
length(orth_dat$sp1[grep("Tgut",orth_dat$sp1)])

length(unique(orth_dat$sp2[grep("Tgut",orth_dat$sp2)]))
length(orth_dat$sp2[grep("Tgut",orth_dat$sp2)])

# Which Ccornix genes are duplicated
# these will represent unique Ccornix genes that have an equally good hit to two 
# or more tgut genes.
ccornix_dupl<-ccornix_sp1[duplicated(ccornix_sp1)]
ccornix_dupl2<-ccornix_sp2[duplicated(ccornix_sp2)]
length(unique(ccornix_dupl))
length(unique(ccornix_dupl2))

# Make non-redundant list of duplicates
dupl<-unique(c(as.character(ccornix_dupl2),
               as.character(ccornix_dupl)))
length(dupl)

head(orth_dat[orth_dat$sp1 %in% dupl,],n = 15)
head(orth_dat[orth_dat$sp2 %in% dupl,])

#
# Remove the duplicated cmon names
orth_dat2<-orth_dat[!(orth_dat$sp1 %in% dupl),]
nrow(orth_dat2)
orth_dat2<-orth_dat2[!(orth_dat2$sp2 %in% dupl),]
nrow(orth_dat2)

# Which tgut genes are duplicated
# these will represent unique Ccornix genes that have an equally good hit to two 
# or more tgut genes.
tgut_sp1<-orth_dat2$sp1[grep("Tgut",orth_dat2$sp1)]
tgut_sp2<-orth_dat2$sp2[grep("Tgut",orth_dat2$sp2)]
length(tgut_sp1)
length(tgut_sp2)
tgut_dupl<-tgut_sp1[duplicated(tgut_sp1)]
tgut_dupl2<-tgut_sp2[duplicated(tgut_sp2)]
length(unique(tgut_dupl))
length(unique(tgut_dupl2))

# Make non-redundant list of duplicates
dupl<-unique(c(as.character(tgut_dupl2),
                 as.character(tgut_dupl)))
length(dupl) # 7 Duplicated genes 

nrow(orth_dat2[orth_dat2$sp1 %in% dupl,])
nrow(orth_dat2[orth_dat2$sp2 %in% dupl,])

# Remove the duplicated tgut names
orth_dat3<-orth_dat2[!(orth_dat2$sp1 %in% dupl),]
nrow(orth_dat3)
orth_dat3<-orth_dat3[!(orth_dat3$sp2 %in% dupl),]
nrow(orth_dat3)


# number of unique 1-1 hits
nrow(orth_dat3)
head(orth_dat3)
orth_dat3[orth_dat3$sp1 == "Tgut|LYN",]
ccornix_3sp1<-orth_dat3[grep("Ccornix",orth_dat3$sp1),]
ccornix_3sp2<-orth_dat3[grep("Ccornix",orth_dat3$sp2),]

tgut_3sp1<-orth_dat3[grep("Tgut",orth_dat3$sp1),]
tgut_3sp2<-orth_dat3[grep("Tgut",orth_dat3$sp2),]

nrow(ccornix_3sp1)
ccornix_3sp1[duplicated(ccornix_3sp1$sp2),]

nrow(ccornix_3sp2)
ccornix_3sp2[duplicated(ccornix_3sp2$sp1),]
# ####

#
# Ortholog groups: After performing the MCL steps
# ####
orth_dat <- read.table("~/RData/ccornix_v_tgut_orthologs_abc.tab",
                       header = FALSE,sep = "\t")
colnames(orth_dat)<-c("orthogroup","genes")
orth_dat$group <- str_pad(seq(1:nrow(orth_dat)),5,pad="0")
orth_dat$orthogroup <- paste(orth_dat$orthogroup,orth_dat$group,":",sep="")
head(orth_dat)

# count:
# the number of genes within each orthogroup
# the number of species within each group
# which species are in each group
orth_dat$n_genes <- vector(length = nrow(orth_dat))
orth_dat$n_sp <- vector(length = nrow(orth_dat))
orth_dat$w_sp <- vector(length = nrow(orth_dat))
for(i in seq(1,nrow(orth_dat),1)){
  genes<-orth_dat$genes[i]
  genes<-as.vector(unlist(str_split(genes,pattern = ",")))
  ngenes<-length(genes)
  nsp<-length(unique(str_replace(genes,pattern = "\\|.*",replacement = "")))
  wsp<-paste(sort(unique(str_replace(genes,pattern = "\\|.*",
                                     replacement = ""))),
             collapse = "-",sep="-")
  orth_dat$n_genes[i] <- ngenes
  orth_dat$n_sp[i] <- nsp
  orth_dat$w_sp[i] <- wsp
}
head(orth_dat)
tail(orth_dat)
nrow(orth_dat)

#plot nr genes per group
orth_hist <- ggplot()+
  geom_histogram(data = orth_dat,aes(n_genes)) +
  xlab("Nr. genes per group")+
  ylab("Count")
orth_hist+my.theme +theme(axis.text = element_text(size = 12))

#plot nr species per group
orth_hist <- ggplot()+
  geom_histogram(data = orth_dat,aes(n_sp)) +
  xlab("Nr. species per group")+
  ylab("Count")
orth_hist+my.theme +theme(axis.text = element_text(size = 12))

# subset only those with all three species
# how many?
nrow(orth_dat[orth_dat$w_sp == "Ccornix-Tgut",])
# how many 1-1 orthologs(i.e. both species, only 2 genes)?
nrow(orth_dat[orth_dat$w_sp == "Ccornix-Tgut" & 
                orth_dat$n_genes == 2,])
one_t_one_orths <- orth_dat[orth_dat$w_sp == "Ccornix-Tgut" & 
                              orth_dat$n_genes == 2,]

tail(one_t_one_orths)
# ####
#
# Subsetting and writing lists
# ####
# Load the common genes
common_5sp <- read.table("~/PhD/crows_project/data/genes/5sp_common_genes.txt",
                         header = FALSE)
common_6sp <- read.table("~/PhD/crows_project/data/genes/6sp_common_genes.txt",
                         header = FALSE)
common_7sp <- read.table("~/PhD/crows_project/data/genes/7sp_common_genes.txt",
                         header = FALSE)
common_8sp <- read.table("~/PhD/crows_project/data/genes/8sp_common_genes.txt",
                         header = FALSE)


# Subset the common genes to include only those that have un-ambiguous orthologs
# Ccornix_3sp1
length(ccornix_3sp2$sp1)
common_5sp_unamb<-common_5sp[
  paste("Ccornix|",common_5sp$V1,sep="") %in% ccornix_3sp2$sp2,]
length(common_5sp_unamb)

common_6sp_unamb<-common_6sp[
  paste("Ccornix|",common_6sp$V1,sep="") %in% ccornix_3sp2$sp2,]
length(common_6sp_unamb)

common_7sp_unamb<-common_7sp[
  paste("Ccornix|",common_7sp$V1,sep="") %in% ccornix_3sp2$sp2,]
length(common_7sp_unamb)

common_8sp_unamb<-common_8sp[
  paste("Ccornix|",common_8sp$V1,sep="") %in% ccornix_3sp2$sp2,]
length(common_8sp_unamb)

# Write "common_*_unamb" to a file
write.table(common_5sp_unamb,
            file = "~/PhD/crows_project/data/genes/5sp_common_unamb_genes.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

write.table(common_6sp_unamb,
            file = "~/PhD/crows_project/data/genes/6sp_common_unamb_genes.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

write.table(common_7sp_unamb,
            file = "~/PhD/crows_project/data/genes/7sp_common_unamb_genes.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

write.table(common_8sp_unamb,
            file = "~/PhD/crows_project/data/genes/8sp_common_unamb_genes.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

# Print a table of which Ccornix gene is which Tgut gene
common_5sp_unamb_ccornix_v_tgut<-ccornix_3sp2[
  ccornix_3sp2$sp2 %in% paste("Ccornix|",common_5sp$V1,sep=""),]
nrow(common_5sp_unamb_ccornix_v_tgut)
head(common_5sp_unamb_ccornix_v_tgut)
write.table(common_5sp_unamb_ccornix_v_tgut[,c(1,2)],
            file = "~/PhD/crows_project/data/genes/5species_ccornix_v_tgut_unamb_genes.tab",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE, sep = "\t")

common_6sp_unamb_ccornix_v_tgut<-ccornix_3sp2[
  ccornix_3sp2$sp2 %in% paste("Ccornix|",common_6sp$V1,sep=""),]
nrow(common_6sp_unamb_ccornix_v_tgut)
write.table(common_6sp_unamb_ccornix_v_tgut[,c(1,2)],
            file = "~/PhD/crows_project/data/genes/6species_ccornix_v_tgut_unamb_genes.tab",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE, sep = "\t")

common_7sp_unamb_ccornix_v_tgut<-ccornix_3sp2[
  ccornix_3sp2$sp2 %in% paste("Ccornix|",common_7sp$V1,sep=""),]
nrow(common_7sp_unamb_ccornix_v_tgut)
write.table(common_7sp_unamb_ccornix_v_tgut[,c(1,2)],
            file = "~/PhD/crows_project/data/genes/7species_ccornix_v_tgut_unamb_genes.tab",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE, sep = "\t")

common_8sp_unamb_ccornix_v_tgut<-ccornix_3sp2[
  ccornix_3sp2$sp2 %in% paste("Ccornix|",common_8sp$V1,sep=""),]
nrow(common_8sp_unamb_ccornix_v_tgut)
write.table(common_8sp_unamb_ccornix_v_tgut[,c(1,2)],
            file = "~/PhD/crows_project/data/genes/8species_ccornix_v_tgut_unamb_genes.tab",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE, sep = "\t")

# ####

#
#
# all-v-all: For CAFE analysis
# Protein sequences Stats ####
Ccornix_gene_dat <- read.table(
  "~/RData/Ccornix_reference_CDS_pep_stats.fa",header = TRUE,
                       sep = "\t")
Ccornix_gene_dat$species<-rep("C. cornix",nrow(Ccornix_gene_dat))
head(Ccornix_gene_dat)

Cmoneduloides_gene_dat <- read.table(
  "~/RData/Cmoneduloides_reference_CDS_pep_stats.fa",header = TRUE,
  sep = "\t")
Cmoneduloides_gene_dat$species<-rep("C. moneduloides",nrow(Cmoneduloides_gene_dat))
head(Cmoneduloides_gene_dat)

Tgut_gene_dat <- read.table(
  "~/RData/Tgut_reference_CDS_pep_stats.fa",header = TRUE,
  sep = "\t")
Tgut_gene_dat$species<-rep("T. guttata",nrow(Tgut_gene_dat))
head(Tgut_gene_dat)

dat<-rbind(Ccornix_gene_dat,Cmoneduloides_gene_dat,Tgut_gene_dat)

prot_l_hist<-ggplot()+
  geom_histogram(data=dat,aes(length,fill=species),position = "dodge")+
  facet_grid(species~.)
prot_l_hist

summary(Ccornix_gene_dat)
summary(Cmoneduloides_gene_dat)
summary(Tgut_gene_dat)
# ####

#
# Analysis of orthAgogue + mcl output: All species####
orth_dat <- read.table("RData/crows_all_orths_abc.txt",
                       header = FALSE,sep = "\t")
colnames(orth_dat)<-c("orthogroup","genes")
orth_dat$group <- str_pad(seq(1:nrow(orth_dat)),5,pad="0")
orth_dat$orthogroup <- paste(orth_dat$orthogroup,orth_dat$group,":",sep="")

orth_dat <- read.table("RData/crows_all_orths_orthofinder.txt",
                       header = FALSE,sep = ":")
colnames(orth_dat)<-c("orthogroup","genes")

head(orth_dat)
tail(orth_dat)
# count:
# the number of genes within each orthogroup
# the number of species within each group
# which species are in each group
orth_dat$n_genes <- vector(length = nrow(orth_dat))
orth_dat$n_sp <- vector(length = nrow(orth_dat))
orth_dat$w_sp <- vector(length = nrow(orth_dat))
for(i in seq(1,nrow(orth_dat),1)){
  genes<-orth_dat$genes[i]
  genes<-as.vector(unlist(str_split(genes,pattern = ",")))
  ngenes<-length(genes)
  nsp<-length(unique(str_replace(genes,pattern = "\\|.*",replacement = "")))
  wsp<-paste(sort(unique(str_replace(genes,pattern = "\\|.*",
                                replacement = ""))),
             collapse = ",",sep=",")
  orth_dat$n_genes[i] <- ngenes
  orth_dat$n_sp[i] <- nsp
  orth_dat$w_sp[i] <- wsp
}
head(orth_dat)
tail(orth_dat)
nrow(orth_dat)

# count: the number of genes for each species
orth_dat$ngenes_p_sp <- vector(length = nrow(orth_dat))
species<-c("Ccornix","Tgut","Cmon")
sp
orth_dat$Ccornix<-vector(length = nrow(orth_dat))
orth_dat$Cmon<-vector(length = nrow(orth_dat))
orth_dat$Tgut<-vector(length = nrow(orth_dat))
for(i in seq(1,nrow(orth_dat),1)){
  for(sp in species){ 
    genes<-orth_dat$genes[i]
    genes<-as.vector(unlist(str_split(genes,pattern = ",")))
    n<-length(grep(sp,genes))
    orth_dat[[sp]][i]<-n
  }
}
orth_dat$description<-rep("OrthoFinder Orthogroup",nrow(orth_dat))
head(orth_dat)
names(orth_dat)
orth_dat_sub<-orth_dat[c(10,1,7,8,9)]
colnames(orth_dat_sub)<-c("Description","ID","Ccornix","Cmon","Tgut")
head(orth_dat_sub)
tail(orth_dat_sub)
orth_dat_sub_m<-melt(orth_dat_sub,id.vars = c("Description","ID"),
                     measure.vars = c("Ccornix","Cmon","Tgut"))
colnames(orth_dat_sub_m)<-c("Description","ID","Species","NGenes")
head(orth_dat_sub_m)

#plot nr genes per group
orth_hist <- ggplot()+
  geom_histogram(data = orth_dat,aes(n_genes)) +
  xlab("Nr. genes per group")+
  ylab("Count")
orth_hist+my.theme +theme(axis.text = element_text(size = 12))

#plot nr species per group
orth_hist <- ggplot()+
  geom_histogram(data = orth_dat,aes(n_sp)) +
  scale_x_continuous(breaks=seq(0,3,1))+
  xlab("Nr. species per group")+
  ylab("Count")
orth_hist+my.theme +theme(axis.text = element_text(size = 12))

# Plot the number of genes per group for each species
orth_hist <- ggplot()+
  geom_bar(data = orth_dat_sub_m,aes(ID,NGenes),stat="identity") +
  xlab("OrthoGroup")+
  ylab("Nr. Genes")+
  facet_grid(Species~.)
orth_hist+my.theme +theme(axis.text = element_text(size = 10),
                          axis.text.x = element_blank(),
                          strip.text.y=element_text(size = 10,angle =0),
                          axis.title = element_text(size = 12),
                          axis.ticks.x = element_blank(),
                          panel.grid = element_blank())


# Filter the table
# Remove gene families that only have members in 1 species
head(orth_dat_sub)
orth_dat_sub_filt<-data.frame(Description = vector(),
                              ID=vector(),
                              Ccornix=vector(),
                              Cmon=vector(),
                              Tgut=vector())
for(fam in orth_dat_sub$ID){
  counts <-as.numeric(orth_dat_sub[
    orth_dat_sub$ID == fam,c(3,4,5)])
  if(length(counts[counts!=0]) == 3){
    #This is a family with members from at all species
    orth_dat_sub_filt<-rbind(orth_dat_sub_filt,orth_dat_sub[orth_dat_sub$ID == fam,])    
  }
}
head(orth_dat_sub_filt)
nrow(orth_dat_sub_filt)
# output the table for CAFE
orth_dat_sub$ID<-gsub(":","",orth_dat_sub$ID)

write.table(orth_dat_sub_filt,
            file = "~/PhD/crows_project/data/genes/OrthoFinder_lCDS_all_crow_species_orthogroup_counts_filt.txt",
            sep ="\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
# ####
# Load the p-values from CAFE.
head(orth_dat_sub,n = 20)



#
# Ccornix v. Cmon
# Checking overlap in Cmon liftover coordinates and Ccornix annotations ####
# Read in Cmon liftover coordinates
cmon <- read.table("RData/Cmon2CcornixV2_CDS_locs.tab", header = FALSE)
colnames(cmon)<-c("chr","start","end","gene")
head(cmon)
# Collapse the CDS to get full gene coordinates
cmon_red<-ddply(cmon,"gene",summarise,
                start =min(start),
                end=max(end),
                chr=unique(chr)[1],
                nchr = length(unique(chr)))
head(cmon_red)
nrow(cmon_red)
length(unique(cmon_red$gene))
cmon[cmon$gene == "augustus_masked-scaffold_0-processed-gene-16.13",]

# Load Ccornix CDS locations
ccornix <- read.table("RData/Ccornix_nomito_lCDS_cdsonly.tab", header = FALSE)
colnames(ccornix)<-c("chr","start","end","gene")

#gffreadCDS <- read.table("~/PhD/crows_project/data/annotation/Ccornix_gffread_CDSs.tab", header = FALSE,sep = "\t")
#colnames(gffreadCDS)<-c("sp","rna","gene")
#head(gffreadCDS)
#nrow(gffreadCDS)
#nrow(ccornix)
#length(unique(gffreadCDS$gene))
#V1<-as.character(unique(ccornix$gene))
#V2<-as.character(gffreadCDS$gene)
#notin<-Reduce(setdiff,list(V1,V2))
#length(notin)

# Collapse the CDS to get full gene coordinates
ccornix_red<-ddply(ccornix,.(gene),summarise,
                start=min(start),
                end=max(end),
                chr=unique(chr)[1],
                nchr = length(unique(chr)))
nrow(ccornix_red)
length(unique(ccornix_red$gene))

# How many Cmon genes have start coords < 0:
nrow(cmon_red[cmon_red$start < 0,])
head(cmon_red[cmon_red$start < 0,])

# Remove genes with start coords < 0
cmon_red2<-cmon_red[cmon_red$start >= 0,]

# Compare the two coordinate lists for overlaps
head(cmon_red[cmon_red$chr == "scaffold_0" 
              & cmon_red$start >= 35497
              & cmon_red$end <= 63507,])

cmon_ccornix_ovlap <- data.frame(
  gene=vector(),start=vector(),end=vector(),chr=vector(),
  cmon_gene=vector(),cmon_start=vector(),cmon_end=vector()
)

# Use start -10 and end + 10 to be conservative.
for(gene in cmon_red$gene){
  start <- cmon_red$start[cmon_red$gene == gene]
  end <- cmon_red$end[cmon_red$gene == gene]
  chr <- as.character(cmon_red$chr[cmon_red$gene == gene])
  ovlap_dat <- ccornix_red[ccornix_red$start >= start-10 &
                             ccornix_red$end <= end+10 &
                             ccornix_red$chr == chr,]
    if(nrow(ovlap_dat) > 0){
      ovlap_dat$cmon_gene <- gene
      ovlap_dat$cmon_start <- start
      ovlap_dat$cmon_end <- end
      #print(ovlap_dat)
      cmon_ccornix_ovlap <- rbind(cmon_ccornix_ovlap,ovlap_dat)
  }
}
cmon_ccornix_ovlap$len <- cmon_ccornix_ovlap$end - cmon_ccornix_ovlap$start
cmon_ccornix_ovlap$cmon_len <- cmon_ccornix_ovlap$cmon_end - cmon_ccornix_ovlap$cmon_start
cmon_ccornix_ovlap$lendiff <- cmon_ccornix_ovlap$len - cmon_ccornix_ovlap$cmon_len
head(cmon_ccornix_ovlap)
# Plot the lengths against each other
lenplot <- ggplot()+
  geom_point(data = cmon_ccornix_ovlap,
             aes(cmon_len/1000000,len/1000000)) +
  xlab("C. moneduloides length (Mb)") +
  ylab("C. cornix length (Mb)")
lenplot

# Plot the lengths against each other
lenplot <- ggplot()+
  geom_point(data = cmon_ccornix_ovlap[
    cmon_ccornix_ovlap$cmon_len < 250000,],
             aes(cmon_len/1000000,len/1000000)) +
  xlab("C. moneduloides length (Mb)") +
  ylab("C. cornix length (Mb)")
lenplot

lenplot <- ggplot()+
  geom_histogram(data = cmon_ccornix_ovlap,
    aes(lendiff/1000000),binwidth = 1)+
  scale_x_continuous(breaks=seq(-50,5,1)) +
  xlab("Length Difference (bp)")+
  ylab("Count")
lenplot

# Nr genes with cmon_start < 0
nrow(cmon_ccornix_ovlap[cmon_ccornix_ovlap$cmon_start < 0,])
head(cmon_ccornix_ovlap[cmon_ccornix_ovlap$cmon_start < 0,])
summary(cmon_ccornix_ovlap[cmon_ccornix_ovlap$cmon_start < 0,])
# those with Cmon start < 0 tend to have a huge length difference with cmon
# lengths being longer than ccornix lengths

# Remove those genes with cmon_start < 0
cmon_ccornix_ovlap<-cmon_ccornix_ovlap[cmon_ccornix_ovlap$cmon_start > 0,]
nrow(cmon_ccornix_ovlap[cmon_ccornix_ovlap$cmon_start < 0,])

# Nr overlaps
nrow(cmon_ccornix_ovlap)
# Nr unique genes in Ccornix that overlap with something in Cmon
length(unique(cmon_ccornix_ovlap$gene))
# Nr unique genes in Cmon that overlap with something in Ccornix
length(cmon_ccornix_ovlap$cmon_gene)
length(unique(cmon_ccornix_ovlap$cmon_gene))

summary(cmon_ccornix_ovlap)

# How many have a |length difference| < 10
nrow(cmon_ccornix_ovlap[abs(cmon_ccornix_ovlap$lendiff) < 15000,])
core_ovlaps <- cmon_ccornix_ovlap[abs(cmon_ccornix_ovlap$lendiff) < 15000,]

lenplot <- ggplot()+
  geom_point(data = core_ovlaps,
    aes(cmon_len/1000000,len/1000000)) +
  xlab("C. moneduloides length (Mb)") +
  ylab("C. cornix length (Mb)")
lenplot

# Any duplications in these genes?
# Nr unique genes in Ccornix that overlap with something in Cmon
length(core_ovlaps$gene)
length(unique(core_ovlaps$gene))
# Which are duplicated
duplicates<-core_ovlaps$gene[duplicated(core_ovlaps$gene)]
duplicates
core_ovlaps[core_ovlaps$gene %in% duplicates,]

# Nr unique genes in Cmon that overlap with something in Ccornix
length(unique(cmon_ccornix_ovlap$cmon_gene[abs(cmon_ccornix_ovlap$lendiff) < 10]))

# Find which Cmon genes have more than one overlap
cmon_duplicates<-unique(cmon_ccornix_ovlap$cmon_gene[duplicated(cmon_ccornix_ovlap$cmon_gene)])
length(cmon_duplicates)
cmon_dupls<-cmon_ccornix_ovlap[cmon_ccornix_ovlap$cmon_gene %in% cmon_duplicates,]

cmon_dupls<-cmon_dupls[order(cmon_dupls$cmon_gene),]
nrow(cmon_dupls)
head(cmon_dupls)
length(unique(cmon_dupls$cmon_gene))

cmon_N_dupls<-cmon_ccornix_ovlap[!(cmon_ccornix_ovlap$cmon_gene %in% cmon_duplicates),]
cmon_N_dupls<-cmon_N_dupls[order(cmon_N_dupls$cmon_gene),]
nrow(cmon_N_dupls)
head(cmon_N_dupls)

# How many Ccornix genes are represented twice in this set
length(unique(cmon_N_dupls$gene)) # some
# Which genes

# Find which ccornix genes have more than one overlap
ccornix_duplicates<-unique(cmon_ccornix_ovlap$gene[duplicated(cmon_ccornix_ovlap$gene)])

ccornix_dupls<-cmon_ccornix_ovlap[cmon_ccornix_ovlap$gene %in% ccornix_duplicates,]
ccornix_dupls<-ccornix_dupls[order(ccornix_dupls$gene),]
nrow(ccornix_dupls)
head(ccornix_dupls)

ccornix_N_dupls<-cmon_ccornix_ovlap[!(cmon_ccornix_ovlap$gene %in% ccornix_duplicates),]
ccornix_N_dupls<-ccornix_N_dupls[order(ccornix_N_dupls$gene),]
nrow(ccornix_N_dupls)
head(ccornix_N_dupls)
# ####

# NicD orthAgogue: PAML- cmon v. tgut ####
nicorths <- read.csv("~/RData/nicd_orthologs_cmon_v_tgut.abc",header = FALSE,
                     sep = "\t")
colnames(nicorths)<-c("sp1","sp2","score")
head(nicorths)
tail(nicorths)
nrow(nicorths)
length(unique(nicorths$sp1))
length(unique(nicorths$sp2))
# Nr Tgut genes
length(unique(nicorths$sp1[grep("Tgutt",nicorths$sp1)]))
length(nicorths$sp1[grep("Tgutt",nicorths$sp1)])
length(unique(nicorths$sp2[grep("Tgutt",nicorths$sp2)]))
length(nicorths$sp2[grep("Tgutt",nicorths$sp2)])

# Nr Cmon genes
length(unique(nicorths$sp1[grep("Cmone",nicorths$sp1)]))
length(nicorths$sp1[grep("Cmone",nicorths$sp1)])

head(nicorths$sp1[grep("Cmone",nicorths$sp1)])
# Subset rows
cmon_sp1<-nicorths$sp1[grep("Cmone",nicorths$sp1)]
cmon_sp2<-nicorths$sp2[grep("Cmone",nicorths$sp2)]
length(cmon_sp1)
length(cmon_sp2)

# Which Cmon genes are duplicated
# these will represent unique Cmon genes that have an equally good hit to two or 
# more tgut genes.
head(cmon_sp1[duplicated(cmon_sp1)])
length(cmon_sp1[duplicated(cmon_sp1)])
cmon_dupl<-cmon_sp1[duplicated(cmon_sp1)]
cmon_dupl2<-cmon_sp2[duplicated(cmon_sp2)]
length(unique(cmon_dupl))
length(unique(cmon_dupl2))

head(cmon_dupl)

head(nicorths[nicorths$sp1 %in% cmon_dupl,])
tail(nicorths[nicorths$sp1 %in% cmon_dupl2,])
# These genes result in nonsense input files to PAML

# Remove the duplicated cmon names
nicorths2<-nicorths[!(nicorths$sp1 %in% cmon_dupl),]

# Which tgut genes are duplicated
# these will represent those Cmon genes that have equally good hits to the 
# same tgut genes.
tgut_sp1<-nicorths2$sp1[grep("Tgutt",nicorths2$sp1)]
tgut_sp2<-nicorths2$sp2[grep("Tgutt",nicorths2$sp2)]
length(tgut_sp1)
length(tgut_sp2)

head(tgut_sp1[duplicated(tgut_sp1)])
length(tgut_sp1[duplicated(tgut_sp1)])
tgut_dupl<-tgut_sp1[duplicated(tgut_sp1)]
tgut_dupl2<-tgut_sp2[duplicated(tgut_sp2)]
length(unique(tgut_dupl))
length(unique(tgut_dupl2))

head(nicorths2[nicorths2$sp1 %in% tgut_dupl,])
tail(nicorths2[nicorths2$sp1 %in% tgut_dupl2,])

# Remove the duplicated tgut names
nicorths2<-nicorths2[!(nicorths2$sp1 %in% tgut_dupl),]

# load Tgut trans_gene tab
trans_gene<-read.table("RData/nicd_tgut_trans_gene.tab",header = TRUE,
                       sep = "\t")

head(trans_gene)
#how many sequences
nrow(trans_gene)
#how many unique CDSs
length(unique(trans_gene$CDS))
# how many unique genes
length(unique(trans_gene$GENE))
# which genes are duplicated
tgut_g_dupl<-unique(trans_gene$GENE[duplicated(trans_gene$GENE)])
length(tgut_g_dupl)
#subset the duplicate data
head(trans_gene[trans_gene$GENE %in% tgut_g_dupl,])
duplicate_cdss<-trans_gene[trans_gene$GENE %in% tgut_g_dupl,]
nrow(duplicate_cdss)
#subset the unique CDSs
unique_cdss<-trans_gene[!(trans_gene$GENE %in% tgut_g_dupl),]
nrow(unique_cdss)

# how many of the duplicates are in the orthologs file
head(paste("Tgutt|",duplicate_cdss$CDS,sep = ""))
head(nicorths2)
head(nicorths2[nicorths2$sp1 %in% paste("Tgutt|",duplicate_cdss$CDS,sep=""),])
nrow(nicorths2[nicorths2$sp1 %in% paste("Tgutt|",duplicate_cdss$CDS,sep=""),])
head(duplicate_cdss[duplicate_cdss$CDS == "ENSTGUT00000011177.1",])
head(duplicate_cdss[duplicate_cdss$GENE == "ENSTGUG00000010709.1",])
# ####





