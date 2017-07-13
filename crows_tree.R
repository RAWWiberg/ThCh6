###
# Plotting and modifying the phylogenetic tree
# Last Modified: 29.04.2017
###
# Clear all
rm(list = ls(all=TRUE))
#install.packages("ade4")
library(ggplot2)
source("~/RData/RScripts/ggplot_theme.R")

library(ape)
library(ade4)
library(phytools)
#source("https://bioconductor.org/biocLite.R")
##?BiocUpgrade
#biocLite()
#biocLite("BiocUpgrade")
#biocLite("ggtree")
library(ggtree)
or <-"#E69F00"
bl <-"#0072B2"

# Load the tree
crow_tree<-read.tree("~/PhD/crows_project/results/phylogeny/RAxML_bestTree.crows_tree_GTRGAMMA_rc40_boot10k")

is.ultrametric(crow_tree)
str(crow_tree)

plot(as.phylo(crow_tree),use.edge.length=FALSE,show.node.label=TRUE)
plot(as.phylo(compute.brlen(crow_tree)))
add.scale.bar()

# Make nice labels for the tree
species_names<-c('paste(bolditalic("C. monedula"))',
                 'paste(bolditalic("C. dauuricus"))',
                 'paste(bolditalic("C. corax"))',
                 'paste(bolditalic("C. hawaiiensis"))',
                 'paste(bolditalic("C. kubaryi"))',
                 'paste(bolditalic("C. frugilegus"))',
                 'paste(bolditalic("C. splendens"))',
                 'paste(bolditalic("C. tasmanicus"))',
                 'paste(bolditalic("C. moneduloides"))',
                 'paste(bolditalic("C. woodfordi"))',
                 'paste(bolditalic("C. brachyrhynchos"))',
                 'paste(bolditalic("C. corone "),bold(Group))',
                 'paste(bolditalic("T. guttata"))')
crow_tree$tip.label<-species_names
tree_data<-data.frame(species=as.factor(crow_tree$tip.label),
                      tool=as.factor(c("n","n","n","y","n","n",
                                       "n","n","y","n","n","n","n")))
pamldata<-data.frame(sp5=as.factor(c("0","1","0","0","0","1","1",
                                          "0","1","0","0","1","1")),
                     sp7=as.factor(c("0","1","1","0","0","1","1",
                                          "1","1","0","0","1","1")),
                     sp8=as.factor(c("0","1","1","1","0","1","1",
                                          "1","1","0","0","1","1")))
rownames(pamldata)<-crow_tree$tip.label

treeplot<- ggtree(crow_tree) %<+% tree_data + 
  geom_tiplab(size=4,align=TRUE,aes(color=tool),
              linesize=0.5,parse=TRUE)+
  scale_colour_manual("",values=c("grey50",or),guide=FALSE)
treeplot<-treeplot+geom_treescale(x = 0.037,offset = -0.5)
gheatmap(treeplot,pamldata,width=0.1,offset = 0.01,
         font.size = 3,colnames_angle=-45,hjust=0)+
  scale_fill_manual(breaks=c("1","0"),values=c("white","black"),guide=FALSE)

# For PAML:
# Load tree again
crow_tree<-read.tree("~/PhD/crows_project/results/phylogeny/RAxML_bestTree.crows_tree_GTRGAMMA_rc40_boot10k")
crow_tree$tip.label
plot(crow_tree)

# Prune the tree
# 5 species set
species5 <- c("3sp_Ccornix","Cdau","Cfru","Csple",
              "Cmon","Tgut")
drop_species5 <- crow_tree$tip.label[which(!(crow_tree$tip.label %in% species5))]
species5_tree <- drop.tip(crow_tree,drop_species5)
plot(as.phylo(unroot(species5_tree)))
write.tree(as.phylo(unroot(species5_tree)),
            "~/PhD/crows_project/results/phylogeny/5species_crows_Nulltree_ur.nex")
# 7 species set
species7 <- c("3sp_Ccornix","Cdau","Cfru","Csple",
              "Cmon","Ctas","Ccorx","Tgut")
drop_species7 <- crow_tree$tip.label[which(!(crow_tree$tip.label %in% species7))]
species7_tree <- drop.tip(crow_tree,drop_species7)
plot(as.phylo(unroot(species7_tree)))
write.tree(as.phylo(unroot(species7_tree)),
           "~/PhD/crows_project/results/phylogeny/7species_crows_Nulltree_ur.nex")
# 8 species set
species8 <- c("3sp_Ccornix","Cdau","Cfru","Csple",
              "Cmon","Ctas","Ccorx","Chaw","Tgut")
drop_species8 <- crow_tree$tip.label[which(!(crow_tree$tip.label %in% species8))]
species8_tree <- drop.tip(crow_tree,drop_species8)
plot(as.phylo(unroot(species8_tree)))
write.tree(as.phylo(unroot(species8_tree)),
           "~/PhD/crows_project/results/phylogeny/8species_crows_Nulltree_ur.nex")


# For CAFE:
# Make dataset of dates for some nodes (C. cornix v C. moneduloides node)
dates<-data.frame(node=c(14,19),age.min=c(36,10),age.max=c(50,11),soft.bounds=c(NA))
# Make the tree ultrametric
crow_tree_um<-chronos(crow_tree,lambda = 0.5,calibration = dates)
plot(crow_tree_um)
add.scale.bar()
crow_tree_um$edge.length<-round(crow_tree_um$edge.length,0)
tips<-c("Tgut","3sp_Ccornix","Cmon")
crow_tree_um <- drop.tip(crow_tree_um,
                         tip = crow_tree_um$tip.label[
                           which(!(crow_tree_um$tip.label %in% tips))])

str(crow_tree_um)
ggtree(crow_tree_um)+scale_x_ggtree()
add.scale.bar(x = 1,y=1.2,lwd = 2)
text("(Million Years)",y=1,x=4,cex = 0.6)

str(crow_tree_um)
# Write the ultrametric pruned tree to a file.
write.tree(crow_tree_um,
           "~/PhD/crows_project/results/darren_phylogeny/crows_CAFEtree.new")


