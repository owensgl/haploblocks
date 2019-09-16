library(tidyverse)
library(ggrepel)

directory <- "/media/owens/Copper/wild_gwas/annuus"

inversion_regions <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt")
chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")
window_size <- 2000000


ann_inversions <- chr_lengths %>% 
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
  filter(species == "Annuus")


ann.fst <- read_tsv(paste0(directory,"/Annuus.tranche90.snp.env.90.bi.remappedHa412HO.TEXvCEN.fst.txt.gz"))



ann.fst.cum <- chr_lengths %>% 
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(ann.fst, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( poscum=pos+tot)

axisdf = ann.fst.cum %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )

ann.fst.cum %>%
  filter(N1 > 5 & N2 > 5) %>%
  filter(Fst > 0.05) %>%
  ggplot(.) +
  
  # Show all points
  geom_point( aes(x=poscum, y=Fst,color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  #Draw mean line:
  geom_line(data=ann.fst.cum %>% 
              mutate(window = floor(poscum/window_size),window_middle=(window*window_size)+(window_size/2)) %>%
              group_by(chr,window,window_middle) %>%
              summarize(meanFst = sum(FstNum)/sum(FstDenom)),
            aes(x=window_middle,y=meanFst),color="black") +
  geom_rect(data=ann_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0.1,ymax=0)) +
  
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


ann.inv.fst <- read_tsv(paste0(directory,"/Ha412HO_inv.v3.pcasites.TEXvCEN.fst.txt"))
ann.noinv.fst <- read_tsv(paste0(directory,"/Annuus.tranche90.snp.env.90.bi.remappedHa412HO.v3noinv.ldfilter.TEXvCEN.fst.txt.gz"))
ann.inv.fst$percentile <- "NA"
for (i in 1:nrow(ann.inv.fst)){
  n_fst <- ann.inv.fst$Fst[i]
  
  percentile <- ann.noinv.fst %>%
    filter(N1 > 5 & N2 > 5) %>%
    mutate(place = case_when(Fst >= n_fst ~ "Higher",
                             TRUE ~ "Lower")) %>%
    group_by(place) %>%
    count() %>%
    ungroup() %>%
    mutate(freq = n / sum(n)) %>%
    filter(place == "Higher") %>%
    pull(freq)
  ann.inv.fst$percentile[i] <- as.numeric(percentile)
}

ann.inv.fst$percent_label <- round(as.numeric(ann.inv.fst$percentile)*100,4)
ann.noinv.fst %>%
  filter(N1 > 5 & N2 > 5) %>%
  ggplot(.,aes(Fst)) + geom_density(fill="grey") +
  geom_point(data=ann.inv.fst,aes(x=Fst,y=0)) +
  geom_text_repel(data=ann.inv.fst %>% filter(percent_label < 5),aes(x=Fst,y=0,label=paste0(gsub("Ha412HO","", chr),"\n",percent_label,"%")),
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
  ggtitle("Annuus Texanus vs Central Fst")



####Both together
pdf("texanus_fst_v2.pdf",height=5,width=12)
ann.fst.cum %>%
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
  geom_line(data=ann.fst.cum %>% 
              filter(N1 > 5 & N2 > 5) %>%
              mutate(window = floor(poscum/window_size),window_middle=(window*window_size)+(window_size/2)) %>%
              group_by(chr,window,window_middle) %>%
              summarize(meanFst = sum(FstNum)/sum(FstDenom),n_sites = n()) %>%
              filter(n_sites > 10),
            aes(x=window_middle,y=meanFst),color="black",alpha=0.8) +
  geom_rect(data=ann_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0.02,ymax=0),fill="black") +
  
  scale_fill_manual(values = rep(c("grey", "skyblue"), 22 )) +
  xlab("Chromosome") + ylab("Fst") +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())

ann.noinv.fst %>%
  filter(N1 > 5 & N2 > 5) %>%
  ggplot(.,aes(Fst)) + geom_density(fill="grey") +
  geom_point(data=ann.inv.fst,aes(x=Fst,y=0)) +
  geom_text_repel(data=ann.inv.fst %>% filter(percent_label < 5),aes(x=Fst,y=0,label=paste0(gsub("Ha412HO","", chr),"\n",percent_label,"%")),
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
  ggtitle("H. annuus vs H. annuus texanus Fst")

dev.off()


