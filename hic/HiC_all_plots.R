library(data.table)
library(tidyverse)
source('findInversions.R')
library(grid)
library(gridExtra)

hic_col_names <- c("chr2","bin2","chr1","bin1","nlinks","genotype")


data_ann <- fread('zcat /media/owens/Copper/wild_gwas/HiC/homer_adapttrimmed/HI.4596.002.HA412.1MB.simple.intertable.txt.gz')
colnames(data_ann) <- hic_col_names
data_ann[,genotype:=NULL]
colnames(data_ann) <- c("chr2","bin2","chr1","bin1","ann_links")


data_pet_dune <- fread('zcat /media/owens/Copper/wild_gwas/HiC/homer_adapttrimmed/HI.5114.004.Illumina_Indexing-4.DTG_HiC_1000_PET_dune.1MB.simple.intertable.txt.gz')
colnames(data_pet_dune) <- hic_col_names
data_pet_nondune <- fread('zcat /media/owens/Copper/wild_gwas/HiC/homer_adapttrimmed/HI.5114.005.Illumina_Indexing-5.DTG_HiC_999_PET_nondune.1MB.simple.intertable.txt.gz')
colnames(data_pet_nondune) <- hic_col_names

data_pet_dune[,genotype:=NULL]
data_pet_nondune[,genotype:=NULL]
colnames(data_pet_dune) <- c("chr2","bin2","chr1","bin1","dune_links")
colnames(data_pet_nondune) <- c("chr2","bin2","chr1","bin1","nondune_links")
data_pet_dune[,nondune_links:=data_pet_nondune[,nondune_links]]
data_pet_dune[,ann_links:=data_ann[,ann_links]]
data_pet_dune[,link_comparison:=dune_links - nondune_links]
data_pet_dune[,dune_v_ann:=dune_links - ann_links]
data_pet_dune[,nondune_v_ann:=nondune_links - ann_links]
data_pet <- data_pet_dune


data_arg_inv <- fread('zcat /media/owens/Copper/wild_gwas/HiC/homer_adapttrimmed/HI.5114.006.Illumina_Indexing-6.DTG_HiC_998_ARG_inv.1MB.simple.intertable.txt.gz')
colnames(data_arg_inv) <- hic_col_names
data_arg_noninv <- fread('zcat /media/owens/Copper/wild_gwas/HiC/homer_adapttrimmed/HI.5114.007.Illumina_Indexing-3.DTG_HiC_995_ARG_notinv.1MB.simple.intertable.txt.gz')
colnames(data_arg_noninv) <- hic_col_names

data_arg_inv[,genotype:=NULL]
data_arg_noninv[,genotype:=NULL]
colnames(data_arg_inv) <- c("chr2","bin2","chr1","bin1","inv_links")
colnames(data_arg_noninv) <- c("chr2","bin2","chr1","bin1","noninv_links")
data_arg_inv[,ann_links:=data_ann[,ann_links]]
data_arg_inv[,noninv_links:=data_arg_noninv[,noninv_links]]
data_arg_inv[,link_comparison:=inv_links - noninv_links]
data_arg_inv[,arginv_v_ann:=inv_links - ann_links]
data_arg_inv[,argnoninv_v_ann:=noninv_links - ann_links]
data_arg <- data_arg_inv


read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt") %>%
  filter(spe != "niveus") %>%
  mutate(id = paste0(chr,".",mds)) -> inversion_locations

read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt") %>%
  filter(spe != "niveus") %>%
  mutate(id = paste0("Ha412HOChr",sprintf('%02d',chr),".",direction,mds)) %>%
  select(sv_name,id, spe) %>%
  inner_join(inversion_locations)-> inversion_locations

inversion_locations %>%
  select(sv_name,id) %>% unique() -> inversions

ann_inv_to_print <- c("ann01.01","ann13.01","ann15.01","ann16.02")
ann_plots <- list()
x <- 1
for (i in 1:nrow(inversions)){
  
  window_size <- 1000000
  
  chosen_inversion <- inversions$sv_name[i]
  chosen_chr <-  inversion_locations %>%
    filter(sv_name == chosen_inversion) %>% pull(chr) %>%.[1]
  chosen_start <- inversion_locations %>%
    filter(sv_name == chosen_inversion) %>% pull(start) %>% min()
  chosen_end <- inversion_locations %>%
    filter(sv_name == chosen_inversion) %>% pull(end) %>% max() 
  chosen_middle = ((floor(chosen_start/1000000) + floor(chosen_end/1000000))/2)*1000000
  chosen_height = ((floor(chosen_end/1000000) - floor(chosen_start/1000000))/2)
  view_start <- chosen_start - 20000000
  view_end <- chosen_end + 20000000
  point_size = 100/((view_end - view_start)/1000000)
  
  if (!chosen_inversion %in% ann_inv_to_print){
    next
  }
  if (chosen_inversion == "ann16.02"){
  p1 <- as_tibble(data_pet[chr2 ==chosen_chr & chr1 ==chosen_chr & bin2 >= view_start & bin2 <= view_end &
                             bin1 >= view_start & bin1 <= view_end]) %>%  
    mutate(dune_v_ann = case_when(bin1 == bin2 ~ 99,
                                  dune_v_ann > 0.3 ~ 0.3,
                                  dune_v_ann < -0.3 ~ -0.3,
                                  dune_v_ann == 0 ~ 0.00001,
                                  TRUE ~ dune_v_ann)) %>% 
    mutate(dune_v_ann = na_if(dune_v_ann, dune_v_ann == 99)) %>% 
    filter(bin1 >= bin2) %>%
    mutate(x = (bin1 + bin2)/2,y=(bin1-bin2)/2000000) %>%
    group_by(y,x,dune_v_ann) %>%
    expand(count = seq(1:4)) %>%
    mutate(bin_id = paste0(x,"_",y)) %>%
    ungroup() %>%
    mutate(y = case_when(count == 2 ~ y+0.5,
                         count == 4 ~ y - 0.5,
                         TRUE ~ y),
           x = case_when(count == 1 ~ x -(window_size/2),
                         count == 3 ~ x + (window_size/2),
                         TRUE ~ x)) %>%
    mutate(y = case_when(y < 0 ~ 0,
                         TRUE ~ y)) %>%
    ggplot(.,aes()) + geom_polygon(aes(x=x/1000000,y=y,fill=dune_v_ann,group = bin_id),size=0) +
    scale_fill_gradient2(low = "#377EB8", mid = "white",
                         high = "#E41A1C", midpoint = 0,name="HiC link\ndifference",limits=c(-0.3,0.3)) +
    #scale_fill_distiller(palette = "Spectral",name="HiC link\ndifference",limits=c(-0.3,0.3)) +
    annotate("segment",x=floor(chosen_start/1000000),y=0,
             yend=chosen_height,xend=chosen_middle/1000000,
             alpha=1,color="black",size=0.1) +
    annotate("segment",x=floor(chosen_end/1000000),y=0,
             yend=chosen_height,xend=chosen_middle/1000000,
             alpha=1,color="black",size=0.1) +
    theme_linedraw() + ylab("") + xlab("Mbp") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          panel.border = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="none",
          text = element_text(size=8)) +
    ggtitle(paste0(chosen_inversion, " ANN vs PET"))
    ann_plots[[x]] = p1
    x <- x+1
    next
  
  }
  
  
  
  
  p1 <- as_tibble(data_arg[chr2 ==chosen_chr & chr1 ==chosen_chr & bin2 >= view_start & bin2 <= view_end &
                             bin1 >= view_start & bin1 <= view_end]) %>%  
    mutate(arginv_v_ann = case_when(bin1 == bin2 ~ 99,
                                    arginv_v_ann > 0.3 ~ 0.3,
                                    arginv_v_ann < -0.3 ~ -0.3,
                                    arginv_v_ann == 0 ~ 0.00001,
                                    TRUE ~ arginv_v_ann)) %>% 
    mutate(arginv_v_ann = na_if(arginv_v_ann, arginv_v_ann == 99)) %>% 
    filter(bin1 >= bin2) %>%
    mutate(x = (bin1 + bin2)/2,y=(bin1-bin2)/2000000) %>%
    group_by(y,x,arginv_v_ann) %>%
    expand(count = seq(1:4)) %>%
    mutate(bin_id = paste0(x,"_",y)) %>%
    ungroup() %>%
    mutate(y = case_when(count == 2 ~ y+0.5,
                         count == 4 ~ y - 0.5,
                         TRUE ~ y),
           x = case_when(count == 1 ~ x -(window_size/2),
                         count == 3 ~ x + (window_size/2),
                         TRUE ~ x)) %>%
    mutate(y = case_when(y < 0 ~ 0,
                         TRUE ~ y)) %>%
    ggplot(.,aes()) + geom_polygon(aes(x=x/1000000,y=y,fill=arginv_v_ann,group = bin_id),size=0) +
    scale_fill_gradient2(low = "#377EB8", mid = "white",
                         high = "#E41A1C", midpoint = 0,name="HiC link\ndifference",limits=c(-0.3,0.3)) +
    #scale_fill_distiller(palette = "Spectral",name="HiC link\ndifference",limits=c(-0.3,0.3)) +
    annotate("segment",x=floor(chosen_start/1000000),y=0,
             yend=chosen_height,xend=chosen_middle/1000000,
             alpha=1,color="black",size=0.1) +
    annotate("segment",x=floor(chosen_end/1000000),y=0,
             yend=chosen_height,xend=chosen_middle/1000000,
             alpha=1,color="black",size=0.1) +
    theme_linedraw() + ylab("") + xlab("Mbp") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          panel.border = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="none",
          text = element_text(size=8)) +
    ggtitle(paste0(chosen_inversion, "ANN vs ARG"))
  ann_plots[[x]] = p1
  x <- x+1
  next
}

dev.off()


read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt") %>%
  filter(spe == "argophyllus") %>%
  mutate(id = paste0(chr,".",mds)) -> arg_inversion_locations
arg_inversion_locations %>%
  select(id) %>% unique() -> arg_inversions
read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt") %>%
  filter(spe == "argophyllus") %>%
  mutate(id = paste0("Ha412HOChr",sprintf('%02d',chr),".",direction,mds)) %>%
  select(sv_name,id) %>%
  inner_join(arg_inversions)-> arg_inversions

arg_genotypes <- read_tsv("HiC/arg_inv_genotypes.txt")
#All petiolaris inversions
arg_plots <- list()
x <- 1
for (i in 1:nrow(arg_inversions)){
  window_size <- 1000000
  chosen_inversion <- arg_inversions$id[i]
  mini_id <- arg_inversions$sv_name[i]
  chosen_chr <-  arg_inversion_locations %>%
    filter(id == chosen_inversion) %>% pull(chr) %>%.[1]
  chosen_start <- arg_inversion_locations %>%
    filter(id == chosen_inversion) %>% pull(start) %>% min()
  chosen_start <- chosen_start
  chosen_end <- arg_inversion_locations %>%
    filter(id == chosen_inversion) %>% pull(end) %>% max() 
  chosen_end <- chosen_end
  chosen_middle = ((floor(chosen_start/1000000) + floor(chosen_end/1000000))/2)*1000000
  chosen_height = ((floor(chosen_end/1000000) - floor(chosen_start/1000000))/2)
  view_start <- chosen_start - 20000000
  view_end <- chosen_end + 20000000
  point_size = 100/((view_end - view_start)/1000000)
  inv_genotype <- arg_genotypes %>% filter(id == mini_id) %>%
    mutate(ARG_inv = case_when(ARG_inv == 0 ~ "0/0",
                               ARG_inv == 1 ~ "0/1",
                               TRUE ~ "1/1"))%>%
    pull(ARG_inv)
  
  notinv_genotype <- arg_genotypes %>% filter(id == mini_id) %>% 
    mutate(ARG_notinv = case_when(ARG_notinv == 0 ~ "0/0",
                               ARG_notinv == 1 ~ "0/1",
                               TRUE ~ "1/1"))%>%
    pull(ARG_notinv)
  
  if (mini_id == "arg06.01" | mini_id == "arg10.01"){
    p1 <- as_tibble(data_arg[chr2 ==chosen_chr & chr1 ==chosen_chr & bin2 >= view_start & bin2 <= view_end &
                               bin1 >= view_start & bin1 <= view_end]) %>%  
      mutate(link_comparison = case_when(bin1 == bin2 ~ 99,
                                         link_comparison > 0.3 ~ 0.3,
                                         link_comparison < -0.3 ~ -0.3,
                                         link_comparison == 0 ~ 0.00001,
                                         TRUE ~ link_comparison)) %>% 
      mutate(link_comparison = na_if(link_comparison, link_comparison == 99)) %>% 
      filter(bin1 >= bin2) %>%
      mutate(x = (bin1 + bin2)/2,y=(bin1-bin2)/2000000) %>%
      group_by(y,x,link_comparison) %>%
      expand(count = seq(1:4)) %>%
      mutate(bin_id = paste0(x,"_",y)) %>%
      ungroup() %>%
      mutate(y = case_when(count == 2 ~ y+0.5,
                           count == 4 ~ y - 0.5,
                           TRUE ~ y),
             x = case_when(count == 1 ~ x -(window_size/2),
                           count == 3 ~ x + (window_size/2),
                           TRUE ~ x)) %>%
      mutate(y = case_when(y < 0 ~ 0,
                           TRUE ~ y)) %>%
      ggplot(.,aes()) + geom_polygon(aes(x=x/1000000,y=y,fill=link_comparison,group = bin_id),size=0) +
      scale_fill_gradient2(low = "#377EB8", mid = "white",
                           high = "#E41A1C", midpoint = 0,name="HiC link\ndifference",limits=c(-0.3,0.3)) +
      #scale_fill_distiller(palette = "Spectral",name="HiC link\ndifference",limits=c(-0.3,0.3)) +
      annotate("segment",x=floor(chosen_start/1000000),y=0,
               yend=chosen_height,xend=chosen_middle/1000000,
               alpha=1,color="black",size=0.1) +
      annotate("segment",x=floor(chosen_end/1000000),y=0,
               yend=chosen_height,xend=chosen_middle/1000000,
               alpha=1,color="black",size=0.1) +
      theme_linedraw() + ylab("") + xlab("Mbp") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            panel.border = element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position="none",
            text = element_text(size=8)) +
      ggtitle(paste0(mini_id," ",inv_genotype," - ",notinv_genotype))
    
    plegend <- as_tibble(data_arg[chr2 ==chosen_chr & chr1 ==chosen_chr & bin2 >= view_start & bin2 <= view_end &
                               bin1 >= view_start & bin1 <= view_end]) %>%  
      mutate(link_comparison = case_when(bin1 == bin2 ~ 99,
                                         link_comparison > 0.3 ~ 0.3,
                                         link_comparison < -0.3 ~ -0.3,
                                         link_comparison == 0 ~ 0.00001,
                                         TRUE ~ link_comparison)) %>% 
      mutate(link_comparison = na_if(link_comparison, link_comparison == 99)) %>% 
      filter(bin1 >= bin2) %>%
      mutate(x = (bin1 + bin2)/2,y=(bin1-bin2)/2000000) %>%
      group_by(y,x,link_comparison) %>%
      expand(count = seq(1:4)) %>%
      mutate(bin_id = paste0(x,"_",y)) %>%
      ungroup() %>%
      mutate(y = case_when(count == 2 ~ y+0.5,
                           count == 4 ~ y - 0.5,
                           TRUE ~ y),
             x = case_when(count == 1 ~ x -(window_size/2),
                           count == 3 ~ x + (window_size/2),
                           TRUE ~ x)) %>%
      mutate(y = case_when(y < 0 ~ 0,
                           TRUE ~ y)) %>%
      ggplot(.,aes()) + geom_polygon(aes(x=x/1000000,y=y,fill=link_comparison,group = bin_id),size=0) +
      scale_fill_gradient2(low = "#377EB8", mid = "white",
                           high = "#E41A1C", midpoint = 0,name="HiC link\ndifference",limits=c(-0.3,0.3)) +
      #scale_fill_distiller(palette = "Spectral",name="HiC link\ndifference",limits=c(-0.3,0.3)) +
      annotate("segment",x=floor(chosen_start/1000000),y=0,
               yend=chosen_height,xend=chosen_middle/1000000,
               alpha=1,color="black",size=0.1) +
      annotate("segment",x=floor(chosen_end/1000000),y=0,
               yend=chosen_height,xend=chosen_middle/1000000,
               alpha=1,color="black",size=0.1) +
      theme_linedraw() + ylab("") + xlab("Mbp") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            panel.border = element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      ggtitle(paste0(mini_id," ",inv_genotype," - ",notinv_genotype))
    if (mini_id == "arg06.01"){
      p1 <- p1 +   annotate("segment",x=130,y=0,yend=5,xend=135,color="black",linetype="dotted",size=0.1) +
        annotate("segment",x=140,y=0,yend=5,xend=135,color="black",linetype="dotted",size=0.1)
      p2 <- p2 +   annotate("segment",x=130,y=0,yend=5,xend=135,color="black",linetype="dotted",size=0.1) +
        annotate("segment",x=140,y=0,yend=5,xend=135,color="black",linetype="dotted",size=0.1)
      p3 <- p3 +   annotate("segment",x=130,y=0,yend=5,xend=135,color="black",linetype="dotted",size=0.1) +
        annotate("segment",x=140,y=0,yend=5,xend=135,color="black",linetype="dotted",size=0.1)
    }
    arg_plots[[x]] = p1
    x <- x+1
    next
  }
}




read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt") %>%
  filter(spe == "petiolaris") %>%
  mutate(id = paste0(chr,".",mds)) -> pet_inversion_locations
pet_inversion_locations %>%
  select(id) %>% unique() -> pet_inversions
read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt") %>%
  filter(spe == "petiolaris") %>%
  mutate(id = paste0("Ha412HOChr",sprintf('%02d',chr),".",direction,mds)) %>%
  select(sv_name,id) %>%
  inner_join(pet_inversions)-> pet_inversions

pet_inv_to_print <- c("pet01.01","pet05.01","pet06.01","pet07.01",
                      "pet08.01","pet09.01","pet10.01","pet11.01",
                      "pet14.01","pet16.01","pet16.02","pet17.01",
                      "pet17.03")

pet_genotypes <- read_tsv("HiC/pet_inv_genotypes.txt")
#All petiolaris inversions
pet_plots <- list()
x <- 1

for (i in 1:nrow(pet_inversions)){
  window_size <- 1000000
  chosen_inversion <- pet_inversions$id[i]
  mini_id <- pet_inversions$sv_name[i]
  chosen_chr <-  pet_inversion_locations %>%
    filter(id == chosen_inversion) %>% pull(chr) %>%.[1]
  chosen_start <- pet_inversion_locations %>%
    filter(id == chosen_inversion) %>% pull(start) %>% min()
  chosen_start <- chosen_start
  chosen_end <- pet_inversion_locations %>%
    filter(id == chosen_inversion) %>% pull(end) %>% max() 
  chosen_end <- chosen_end
  chosen_middle = ((floor(chosen_start/1000000) + floor(chosen_end/1000000))/2)*1000000
  chosen_height = ((floor(chosen_end/1000000) - floor(chosen_start/1000000))/2)
  view_start <- chosen_start - 20000000
  view_end <- chosen_end + 20000000
  point_size = 100/((view_end - view_start)/1000000)
  dune_genotype <- pet_genotypes %>% filter(id == mini_id) %>% 
    mutate(PET_dune = case_when(PET_dune == 0 ~ "0/0",
                                PET_dune == 1 ~ "0/1",
                               TRUE ~ "1/1"))%>%
    pull(PET_dune)
  nondune_genotype <- pet_genotypes %>% filter(id == mini_id) %>% 
    mutate(PET_nondune = case_when(PET_nondune == 0 ~ "0/0",
                                   PET_nondune == 1 ~ "0/1",
                               TRUE ~ "1/1"))%>%
    pull(PET_nondune)
  
  if (!mini_id %in% pet_inv_to_print){
    next
  }
  
  
  p1 <- as_tibble(data_pet[chr2 ==chosen_chr & chr1 ==chosen_chr & bin2 >= view_start & bin2 <= view_end &
                             bin1 >= view_start & bin1 <= view_end]) %>%  
    mutate(link_comparison = case_when(bin1 == bin2 ~ 99,
                                       link_comparison > 0.3 ~ 0.3,
                                       link_comparison < -0.3 ~ -0.3,
                                       link_comparison == 0 ~ 0.00001,
                                       TRUE ~ link_comparison)) %>% 
    mutate(link_comparison = na_if(link_comparison, link_comparison == 99)) %>% 
    filter(bin1 >= bin2) %>%
    mutate(x = (bin1 + bin2)/2,y=(bin1-bin2)/2000000) %>%
    group_by(y,x,link_comparison) %>%
    expand(count = seq(1:4)) %>%
    mutate(bin_id = paste0(x,"_",y)) %>%
    ungroup() %>%
    mutate(y = case_when(count == 2 ~ y+0.5,
                         count == 4 ~ y - 0.5,
                         TRUE ~ y),
           x = case_when(count == 1 ~ x -(window_size/2),
                         count == 3 ~ x + (window_size/2),
                         TRUE ~ x)) %>%
    mutate(y = case_when(y < 0 ~ 0,
                         TRUE ~ y)) %>%
    ggplot(.,aes()) + geom_polygon(aes(x=x/1000000,y=y,fill=link_comparison,group = bin_id),size=0) +
    scale_fill_gradient2(low = "#377EB8", mid = "white",
                         high = "#E41A1C", midpoint = 0,name="HiC link\ndifference",limits=c(-0.3,0.3)) +
    #scale_fill_distiller(palette = "RdBu",name="HiC link\ndifference",limits=c(-0.3,0.3)) +
    annotate("segment",x=floor(chosen_start/1000000),y=0,
             yend=chosen_height,xend=chosen_middle/1000000,
             alpha=1,color="black",size=0.1) +
    annotate("segment",x=floor(chosen_end/1000000),y=0,
             yend=chosen_height,xend=chosen_middle/1000000,
             alpha=1,color="black",size=0.1) +
    theme_linedraw() + ylab("") + xlab("Mbp") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          panel.border = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="none",
          text = element_text(size=8)) +
    ggtitle(paste0(mini_id," ",dune_genotype," - ",nondune_genotype))
  
  
  
  
  
  if (mini_id == "pet01.01"){
    p1 <- p1 +   annotate("segment",x=0,y=0,yend=6,xend=6,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=12,y=0,yend=6,xend=6,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=90,y=0,yend=4,xend=94,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=98,y=0,yend=4,xend=94,color="black",linetype="dotted",size=0.1)
    p2 <- p2 +   annotate("segment",x=0,y=0,yend=6,xend=6,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=12,y=0,yend=6,xend=6,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=90,y=0,yend=4,xend=94,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=98,y=0,yend=4,xend=94,color="black",linetype="dotted",size=0.1)
    p3 <- p3 +   annotate("segment",x=0,y=0,yend=6,xend=6,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=12,y=0,yend=6,xend=6,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=90,y=0,yend=4,xend=94,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=98,y=0,yend=4,xend=94,color="black",linetype="dotted",size=0.1)
  }else if (mini_id == "pet09.01"){
    p1 <- p1 +   annotate("segment",x=105,y=0,yend=9,xend=114,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=123,y=0,yend=9,xend=114,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=128,y=0,yend=6,xend=134,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=140,y=0,yend=6,xend=134,color="black",linetype="dotted",size=0.1)
    p2 <- p2 +   annotate("segment",x=105,y=0,yend=9,xend=114,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=123,y=0,yend=9,xend=114,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=128,y=0,yend=6,xend=134,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=140,y=0,yend=6,xend=134,color="black",linetype="dotted",size=0.1)
    p3 <- p3 +   annotate("segment",x=105,y=0,yend=9,xend=114,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=123,y=0,yend=9,xend=114,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=128,y=0,yend=6,xend=134,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=140,y=0,yend=6,xend=134,color="black",linetype="dotted",size=0.1)
  }else if (mini_id == "pet05.01"){
    p1 <- p1 +   annotate("segment",x=154,y=0,yend=13.5,xend=167.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=181,y=0,yend=13.5,xend=167.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=182,y=0,yend=1.5,xend=183.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=185,y=0,yend=1.5,xend=183.5,color="black",linetype="dotted",size=0.1)
    p2 <- p2 +   annotate("segment",x=154,y=0,yend=13.5,xend=167.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=181,y=0,yend=13.5,xend=167.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=182,y=0,yend=1.5,xend=183.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=185,y=0,yend=1.5,xend=183.5,color="black",linetype="dotted",size=0.1)
    p3 <- p3 +   annotate("segment",x=154,y=0,yend=13.5,xend=167.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=181,y=0,yend=13.5,xend=167.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=182,y=0,yend=1.5,xend=183.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=185,y=0,yend=1.5,xend=183.5,color="black",linetype="dotted",size=0.1)
  }else if (mini_id == "pet12.01"){
    p1 <- p1 +   annotate("segment",x=89,y=0,yend=3,xend=92,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=95,y=0,yend=3,xend=92,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=98,y=0,yend=6.5,xend=104.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=111,y=0,yend=6.5,xend=104.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=155,y=0,yend=4.5,xend=159.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=164,y=0,yend=4.5,xend=159.5,color="black",linetype="dotted",size=0.1)
    p2 <- p2 +   annotate("segment",x=89,y=0,yend=3,xend=92,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=95,y=0,yend=3,xend=92,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=98,y=0,yend=6.5,xend=104.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=111,y=0,yend=6.5,xend=104.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=155,y=0,yend=4.5,xend=159.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=164,y=0,yend=4.5,xend=159.5,color="black",linetype="dotted",size=0.1)
    p3 <- p3 +   annotate("segment",x=89,y=0,yend=3,xend=92,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=95,y=0,yend=3,xend=92,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=98,y=0,yend=6.5,xend=104.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=111,y=0,yend=6.5,xend=104.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=155,y=0,yend=4.5,xend=159.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=164,y=0,yend=4.5,xend=159.5,color="black",linetype="dotted",size=0.1)
  }else if (mini_id == "pet16.02"){
    p1 <- p1 +   annotate("segment",x=6,y=0,yend=8.5,xend=14.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=23,y=0,yend=8.5,xend=14.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=34,y=0,yend=5.5,xend=39.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=45,y=0,yend=5.5,xend=39.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=108,y=0,yend=7.5,xend=115.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=123,y=0,yend=7.5,xend=115.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=206,y=0,yend=2,xend=208,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=210,y=0,yend=2,xend=208,color="black",linetype="dotted",size=0.1)   
    p2 <- p2 +   annotate("segment",x=6,y=0,yend=8.5,xend=14.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=23,y=0,yend=8.5,xend=14.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=34,y=0,yend=5.5,xend=39.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=45,y=0,yend=5.5,xend=39.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=108,y=0,yend=7.5,xend=115.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=123,y=0,yend=7.5,xend=115.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=206,y=0,yend=2,xend=208,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=210,y=0,yend=2,xend=208,color="black",linetype="dotted",size=0.1)  
    p3 <- p3 +   annotate("segment",x=6,y=0,yend=8.5,xend=14.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=23,y=0,yend=8.5,xend=14.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=34,y=0,yend=5.5,xend=39.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=45,y=0,yend=5.5,xend=39.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=108,y=0,yend=7.5,xend=115.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=123,y=0,yend=7.5,xend=115.5,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=206,y=0,yend=2,xend=208,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=210,y=0,yend=2,xend=208,color="black",linetype="dotted",size=0.1)  
  }
  pet_plots[[x]] = p1
  x <- x+1
  next
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(plegend)

pdf("HiC/HiC_adaptTrim_NGM_homer10mq_simple.allSV.pdf",height=9.7,width=7.2)
  grid.arrange(ann_plots[[1]],ann_plots[[2]],ann_plots[[3]],ann_plots[[4]],
               arg_plots[[1]],arg_plots[[2]],pet_plots[[1]],pet_plots[[2]],
               pet_plots[[3]],pet_plots[[4]],pet_plots[[5]],pet_plots[[6]],
               pet_plots[[7]],pet_plots[[8]],pet_plots[[9]],pet_plots[[10]],
               pet_plots[[11]],pet_plots[[12]],pet_plots[[13]],mylegend,
               nrow=6)

dev.off()

