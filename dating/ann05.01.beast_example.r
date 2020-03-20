## Plot the tree from BEAST
## Kaichi Huang 2019 Jul

library(ggtree)
library(treeio)
library(ape)
library(tidyverse)
library(RColorBrewer)

mu <- 6.9e-3
my.colors <- c("#000000","#FFCE6E","#338732","#D39117","#4BB7B7","#447499","#645084","#000000")
my.labels <- c(' ','H. annuus type 0','H. argophyllus','H. annuus type 2','H. petiolaris petiolaris','H. petiolaris fallax', 'H. niveus canescens','Perennial outgroup')

tree <- read.beast("Ha412HOChr05.pos1.inv.fixTree.mcc.tre")
#ggtree(tree) + geom_text(aes(label=node))
tree <- groupOTU(tree, list(seq(2,6),c(10,11),c(1,7,8,9,14),c(15,16),c(18,20),c(17,19),c(12,13)))

p <- ggtree(tree, aes(color=group))
revts(p) +
  theme_tree2() +
  geom_range(range="height_0.95_HPD", color="grey", alpha=.6, size=2, branch.length="height") +
  theme(legend.position="right", legend.title=element_blank()) +
  scale_color_manual(values=my.colors, labels=my.labels) +
  scale_x_continuous(limits=c(-0.018, 0.0005), breaks=seq(0,3,0.5)*(-mu), minor_breaks=seq(0,3,0.25)*(-mu), labels=seq(0,3,0.5))
