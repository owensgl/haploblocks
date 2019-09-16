library(tidyverse)
library(ggrepel)

directory <- "/media/owens/Copper/wild_gwas/petiolaris"

inversion_regions <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt")
chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")
window_size <- 2000000


pet_inversions <- chr_lengths %>% 
  select(chr,end) %>%
  rename(chr_len = end) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(inversion_regions, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate( startcum=start+tot,endcum=end+tot) %>%
  filter(species == "Petiolaris")


mon.fst <- read_tsv(paste0(directory,"/Petiolaris.tranche90.snp.petPF.90.bi.remappedHa412HO.MON.DvND.fst.txt.gz"))



mon.fst.cum <- chr_lengths %>% 
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(mon.fst, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( poscum=pos+tot)

axisdf = mon.fst.cum %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )

mon.fst.cum %>%
  filter(N1 > 5 & N2 > 5) %>%
  filter(Fst > 0.05) %>%
  ggplot(.) +
  
  # Show all points
  geom_point( aes(x=poscum, y=Fst,color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  #Draw mean line:
  geom_line(data=mon.fst.cum %>% 
              mutate(window = floor(poscum/window_size),window_middle=(window*window_size)+(window_size/2)) %>%
              group_by(chr,window,window_middle) %>%
              summarize(meanFst = sum(FstNum)/sum(FstDenom)),
            aes(x=window_middle,y=meanFst),color="black") +
  geom_rect(data=pet_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0.1,ymax=0)) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())


mon.inv.fst <- read_tsv(paste0(directory,"/Ha412HO_inv.v3.pcasites.MON.DvND.fst.txt"))
mon.noinv.fst <- read_tsv(paste0(directory,"/Petiolaris.tranche90.snp.petPF.90.bi.remappedHa412HO.v3noinv.ldfilter.MON.DvND.fst.txt.gz"))
mon.inv.fst$percentile <- "NA"
for (i in 1:nrow(mon.inv.fst)){
  n_fst <- mon.inv.fst$Fst[i]
  
  percentile <- mon.noinv.fst %>%
    filter(N1 > 5 & N2 > 5) %>%
    mutate(place = case_when(Fst >= n_fst ~ "Higher",
                             TRUE ~ "Lower")) %>%
    group_by(place) %>%
    count() %>%
    ungroup() %>%
    mutate(freq = n / sum(n)) %>%
    filter(place == "Higher") %>%
    pull(freq)
  mon.inv.fst$percentile[i] <- as.numeric(percentile)
}

mon.inv.fst$percent_label <- round(as.numeric(mon.inv.fst$percentile)*100,4)
mon.noinv.fst %>%
  filter(N1 > 5 & N2 > 5) %>%
  ggplot(.,aes(Fst)) + geom_density(fill="grey") +
  geom_point(data=mon.inv.fst,aes(x=Fst,y=0)) +
  geom_text_repel(data=mon.inv.fst %>% filter(percent_label < 5),aes(x=Fst,y=0,label=paste0(gsub("Ha412HO","", chr),"\n",percent_label,"%")),
                  nudge_y      = 5,
                  direction    = "x",
                  angle        = 0,
                  vjust        = 0,
                  segment.size = 0.2
  ) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  ylab("Density") +
  ggtitle("Texas Dune vs Non-dune Fst")

#######
directory <- "/media/owens/Copper/wild_gwas/petiolaris"
gsd.fst <- read_tsv(paste0(directory,"/Petiolaris.tranche90.snp.petPF.90.bi.remappedHa412HO.GSD.DvND.fst.txt.gz"))

gsd.fst.cum <- chr_lengths %>% 
  select(chr,end) %>%
  rename(chr_len = end) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gsd.fst, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( poscum=pos+tot)

axisdf = gsd.fst.cum %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )

gsd.fst.cum %>%
  filter(N1 > 5 & N2 > 5) %>%
  filter(Fst > 0.05) %>%
  ggplot(.) +
  
  # Show all points
  geom_point( aes(x=poscum, y=Fst,color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  #Draw mean line:
  geom_line(data=gsd.fst.cum %>% 
              mutate(window = floor(poscum/window_size),window_middle=(window*window_size)+(window_size/2)) %>%
              group_by(chr,window,window_middle) %>%
              summarize(meanFst = sum(FstNum)/sum(FstDenom)),
            aes(x=window_middle,y=meanFst),color="black") +
  geom_rect(data=pet_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0.1,ymax=0)) +
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())



gsd.inv.fst <- read_tsv(paste0(directory,"/Ha412HO_inv.v3.pcasites.GSD.DvND.fst.txt"))
gsd.noinv.fst <- read_tsv(paste0(directory,"/Petiolaris.tranche90.snp.petPF.90.bi.remappedHa412HO.v3noinv.ldfilter.GSD.DvND.fst.txt.gz"))
gsd.inv.fst$percentile <- "NA"
for (i in 1:nrow(gsd.inv.fst)){
  n_fst <- gsd.inv.fst$Fst[i]
  
  percentile <- gsd.noinv.fst %>%
    filter(N1 > 5 & N2 > 5) %>%
    mutate(place = case_when(Fst >= n_fst ~ "Higher",
                             TRUE ~ "Lower")) %>%
    group_by(place) %>%
    count() %>%
    ungroup() %>%
    mutate(freq = n / sum(n)) %>%
    filter(place == "Higher") %>%
    pull(freq)
  gsd.inv.fst$percentile[i] <- as.numeric(percentile)
}

gsd.inv.fst$percent_label <- round(as.numeric(gsd.inv.fst$percentile)*100,1)

gsd.noinv.fst %>%
  filter(N1 > 5 & N2 > 5) %>%
  ggplot(.,aes(Fst)) + geom_density(fill="grey") +
  geom_point(data=gsd.inv.fst,aes(x=Fst,y=0)) +
  geom_text_repel(data=gsd.inv.fst %>% filter(percent_label < 5),aes(x=Fst,y=0,label=paste0(gsub("Ha412HO","", chr),"\n",percent_label,"%")),
                  nudge_y      = 5,
                  direction    = "x",
                  angle        = 0,
                  vjust        = 0,
                  segment.size = 0.2
  ) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  ylab("Density") +
  ggtitle("Colorado Dune vs Non-dune Fst")




#######
directory <- "/media/owens/Copper/wild_gwas/petiolaris"
gsdmon.fst <- read_tsv(paste0(directory,"/Petiolaris.tranche90.snp.petPF.90.bi.remappedHa412HO.GSDvMON.fst.txt.gz"))

gsdmon.fst.cum <- chr_lengths %>% 
  select(chr,end) %>%
  rename(chr_len = end) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gsdmon.fst, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( poscum=pos+tot)

axisdf = gsdmon.fst.cum %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )

gsdmon.fst.cum %>%
  filter(Fst > 0.05) %>%
  ggplot(.) +
  
  # Show all points
  geom_point( aes(x=poscum, y=Fst,color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  #Draw mean line:
  geom_line(data=gsdmon.fst.cum %>% 
              mutate(window = floor(poscum/window_size),window_middle=(window*window_size)+(window_size/2)) %>%
              group_by(chr,window,window_middle) %>%
              summarize(meanFst = sum(FstNum)/sum(FstDenom)),
            aes(x=window_middle,y=meanFst),color="black") +
  geom_rect(data=pet_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0.1,ymax=0)) +
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())

  

####Both together
pdf("dune_fst_part1.pdf",height=4,width=8)
gsd.fst.cum %>%
  filter(N1 > 5 & N2 > 5) %>%
  filter(Fst > 0.05) %>%
  ggplot(.) +
  
  # Show all points
  #geom_point( aes(x=poscum, y=Fst,color=as.factor(chr)), alpha=0.8, size=1.3) +
  
  # custom X axis:
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  geom_rect(data=chr_lengths %>% 
              select(chr,start,end) %>%
              rename(chr_len = end) %>%
              # Calculate cumulative position of each chromosome
              mutate(cumstart=cumsum(chr_len)-chr_len, 
                     cumend=cumsum(chr_len)),
            aes(xmin=cumstart,xmax=cumend,ymin=0,ymax=0.75,fill=as.factor(chr)),alpha=0.4) +
  #Draw mean line:
  geom_line(data=gsd.fst.cum %>% 
              filter(N1 > 5 & N2 > 5) %>%
              mutate(window = floor(poscum/window_size),window_middle=(window*window_size)+(window_size/2)) %>%
              group_by(chr,window,window_middle) %>%
              summarize(meanFst = sum(FstNum)/sum(FstDenom),n_sites = n()) %>%
              filter(n_sites > 10),
            aes(x=window_middle,y=meanFst),color="#B7A083",alpha=0.8) +
  geom_line(data=mon.fst.cum %>% 
              filter(N1 > 5 & N2 > 5) %>%
              mutate(window = floor(poscum/window_size),window_middle=(window*window_size)+(window_size/2)) %>%
              group_by(chr,window,window_middle) %>%
              summarize(meanFst = sum(FstNum)/sum(FstDenom),n_sites = n()) %>%
              filter(n_sites > 10),
            aes(x=window_middle,y=meanFst),color="#5C4C39",alpha=0.8) +
  geom_rect(data=pet_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0.02,ymax=0),fill="black") +

  scale_fill_manual(values = rep(c("light grey", "white"), 22 )) +
  xlab("Chromosome") + ylab(expression(F[ST])) +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())
dev.off()
pdf("dune_fst_part2.pdf",height=5,width=5)
mon.noinv.fst %>%
  mutate(type = "Texas") %>%
  rbind(gsd.noinv.fst %>% mutate(type = "Colorado")) %>%
  filter(N1 > 5 & N2 > 5) %>%
  ggplot(.,aes(y=Fst,x=type)) + geom_violin(aes(fill=type)) +
  scale_fill_manual(values=c("#B7A083","#5C4C39")) +
  geom_point(data=mon.inv.fst %>% mutate(type = "Texas"),aes(y=Fst,x=type)) +
  geom_text_repel(data=mon.inv.fst %>% mutate(type = "Texas") %>% filter(percent_label < 5),aes(y=Fst,x="Texas",label=paste0(gsub("Ha412HO","", chr),"\n",percent_label,"%")),
                  nudge_x      = -0.25,
                  direction    = "y",
                  angle        = 0,
                  vjust        = 1,
                  segment.size = 0.2
  ) +
  geom_text_repel(data=gsd.inv.fst %>% mutate(type = "Colorado") %>% filter(percent_label < 5),aes(y=Fst,x="Colorado",label=paste0(gsub("Ha412HO","", chr),"\n",percent_label,"%")),
                  nudge_x      = 0.25,
                  direction    = "y",
                  angle        = 0,
                  vjust        = 1,
                  segment.size = 0.2
  ) + geom_point(data=gsd.inv.fst %>% mutate(type = "Colorado"),aes(y=Fst,x=type)) +
  theme_bw() +
  theme( legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  ylab(expression(paste("Dune Non-dune ", F[S][T],sep=""))) +
  xlab("Location")
dev.off()




