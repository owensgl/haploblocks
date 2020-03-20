## Plot phylogeny illustration using ggtree
## Kaichi Huang 2019 Jul

library(ggtree)
library(treeio)
library(ape)
library(tidyverse)
library(RColorBrewer)

tree <- read.nexus("phylo_130140.contree.nexus") # nexus tree edited and exported by Figtree from .contree

pdf("phylo_130140_example.pdf", width=8, height=10)
tree2 <- groupClade(tree, c(389,431))
ggtree(tree2, aes(color=group)) + 
  scale_color_manual(values=c("#D39117","#0D6D16","#6DBF76"), labels=c('H. annuus','H. argophyllus type 2','H. argophyllus type 0')) +
  theme(legend.position="right", legend.title=element_blank()) + geom_treescale(x=0.1,y=-3,width=0.1,offset=-4,fontsize=3)
dev.off()
