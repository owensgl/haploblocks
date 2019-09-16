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
pdf("HiC/HiC_adaptTrim_NGM_homer10mq_simple.ann.pdf",height=6,width=8)

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
          axis.ticks.y=element_blank()) +
    ggtitle(paste0(chosen_inversion, "ANN vs PET_dune"))
  
  p2 <- as_tibble(data_pet[chr2 ==chosen_chr & chr1 ==chosen_chr & bin2 >= view_start & bin2 <= view_end &
                       bin1 >= view_start & bin1 <= view_end]) %>%  
    mutate(nondune_v_ann = case_when(bin1 == bin2 ~ 99,
                                     nondune_v_ann > 0.3 ~ 0.3,
                                     nondune_v_ann < -0.3 ~ -0.3,
                                     nondune_v_ann == 0 ~ 0.00001,
                                  TRUE ~ nondune_v_ann)) %>% 
    mutate(nondune_v_ann = na_if(nondune_v_ann, nondune_v_ann == 99)) %>% 
    filter(bin1 >= bin2) %>%
    mutate(x = (bin1 + bin2)/2,y=(bin1-bin2)/2000000) %>%
    group_by(y,x,nondune_v_ann) %>%
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
    ggplot(.,aes()) + geom_polygon(aes(x=x/1000000,y=y,fill=nondune_v_ann,group = bin_id),size=0) +
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
    ggtitle(paste0(chosen_inversion, "ANN vs PET_nondune"))
  
  
  p3 <- as_tibble(data_arg[chr2 ==chosen_chr & chr1 ==chosen_chr & bin2 >= view_start & bin2 <= view_end &
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
          axis.ticks.y=element_blank()) +
    ggtitle(paste0(chosen_inversion, "ANN vs ARG_early"))
  
  p4 <- as_tibble(data_arg[chr2 ==chosen_chr & chr1 ==chosen_chr & bin2 >= view_start & bin2 <= view_end &
                       bin1 >= view_start & bin1 <= view_end]) %>%  
    mutate(argnoninv_v_ann = case_when(bin1 == bin2 ~ 99,
                                       argnoninv_v_ann > 0.3 ~ 0.3,
                                       argnoninv_v_ann < -0.3 ~ -0.3,
                                       argnoninv_v_ann == 0 ~ 0.00001,
                                    TRUE ~ argnoninv_v_ann)) %>% 
    mutate(argnoninv_v_ann = na_if(argnoninv_v_ann, argnoninv_v_ann == 99)) %>% 
    filter(bin1 >= bin2) %>%
    mutate(x = (bin1 + bin2)/2,y=(bin1-bin2)/2000000) %>%
    group_by(y,x,argnoninv_v_ann) %>%
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
    ggplot(.,aes()) + geom_polygon(aes(x=x/1000000,y=y,fill=argnoninv_v_ann,group = bin_id),size=0) +
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
    ggtitle(paste0(chosen_inversion, "ANN vs ARG_late"))
  

  



  
  if (chosen_inversion == "arg06.01"){
    p1 <- p1 +   annotate("segment",x=130,y=0,yend=5,xend=135,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=140,y=0,yend=5,xend=135,color="black",linetype="dotted",size=0.1)
    p2 <- p2 +   annotate("segment",x=130,y=0,yend=5,xend=135,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=140,y=0,yend=5,xend=135,color="black",linetype="dotted",size=0.1)
    p3 <- p3 +   annotate("segment",x=130,y=0,yend=5,xend=135,color="black",linetype="dotted",size=0.1) +
      annotate("segment",x=140,y=0,yend=5,xend=135,color="black",linetype="dotted",size=0.1)
    
    }else if (chosen_inversion == "pet01.01"){
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
    }else if (chosen_inversion == "pet09.01"){
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
    }else if (chosen_inversion == "pet05.01"){
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
    }else if (chosen_inversion == "pet12.01"){
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
    }else if (chosen_inversion == "pet16.02"){
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
  print(
    grid.arrange(p1,p2,p3,p4,nrow=2)
  )
}
dev.off()

pdf("HiC/HiC_adaptTrim_NGM_homer10mq_simple.arg6.15.pdf",height=7,width=7)
as_tibble(data_arg[chr2 == "Ha412HOChr06" & chr1 == "Ha412HOChr15"]) %>%  
  ggplot(.,aes()) + geom_tile(aes(x=bin1/1000000,y=bin2/1000000,fill=inv_links)) +
  scale_fill_distiller(palette = "Spectral",name="HiC links") +
  theme_bw() + ylab("Chr6 Mbp") + xlab("Chr15 Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  ggtitle("Arg_inv")

as_tibble(data_arg[chr2 == "Ha412HOChr06" & chr1 == "Ha412HOChr15"]) %>%  
  ggplot(.,aes()) + geom_tile(aes(x=bin1/1000000,y=bin2/1000000,fill=noninv_links)) +
  scale_fill_distiller(palette = "Spectral",name="HiC links") +
  theme_bw() + ylab("Chr6 Mbp") + xlab("Chr15 Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  ggtitle("Arg_noninv")
dev.off()



png("HiC/HiC_adaptTrim_NGM_homer10mq_simple.arg_inv.all.png",height=4000,width=4000)
as_tibble(data_arg[chr2 != chr1]) %>% 
  mutate(chr1 = gsub("Ha412HOChr","",chr1),
         chr2 = gsub("Ha412HOChr","",chr2)) %>%
  filter(chr1 > chr2) %>%
  mutate(inv_links = case_when(inv_links > 0.5 ~ 0.5,
                               TRUE ~ inv_links)) %>%
  ggplot(.,aes()) + geom_tile(aes(x=bin1/1000000,y=bin2/1000000,fill=inv_links)) +
  scale_fill_distiller(palette = "Spectral",name="HiC links") +
  theme_bw() + ylab("Chr6 Mbp") + xlab("Chr15 Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  ggtitle("Arg_inv") +
  facet_grid(chr2~chr1,scales="free")
dev.off()

png("HiC/HiC_adaptTrim_NGM_homer10mq_simple.arg_noninv.all.png",height=4000,width=4000)
as_tibble(data_arg[chr2 != chr1]) %>% 
  mutate(chr1 = gsub("Ha412HOChr","",chr1),
         chr2 = gsub("Ha412HOChr","",chr2)) %>%
  filter(chr1 > chr2) %>%
  mutate(noninv_links = case_when(noninv_links > 0.5 ~ 0.5,
                               TRUE ~ noninv_links)) %>%
  ggplot(.,aes()) + geom_tile(aes(x=bin1/1000000,y=bin2/1000000,fill=noninv_links)) +
  scale_fill_distiller(palette = "Spectral",name="HiC links") +
  theme_bw() + ylab("Chr6 Mbp") + xlab("Chr15 Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  ggtitle("Arg_noninv") +
  facet_grid(chr2~chr1,scales="free")
dev.off()

png("HiC/HiC_adaptTrim_NGM_homer10mq_simple.pet_dune.all.png",height=4000,width=4000)
as_tibble(data_pet[chr2 != chr1]) %>% 
  mutate(chr1 = gsub("Ha412HOChr","",chr1),
         chr2 = gsub("Ha412HOChr","",chr2)) %>%
  filter(chr1 > chr2) %>%
  mutate(dune_links = case_when(dune_links > 0.5 ~ 0.5,
                               TRUE ~ dune_links)) %>%
  ggplot(.,aes()) + geom_tile(aes(x=bin1/1000000,y=bin2/1000000,fill=dune_links)) +
  scale_fill_distiller(palette = "Spectral",name="HiC links") +
  theme_bw() + ylab("Chr6 Mbp") + xlab("Chr15 Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  ggtitle("Pet_dune") +
  facet_grid(chr2~chr1,scales="free")
dev.off()
  
png("HiC/HiC_adaptTrim_NGM_homer10mq_simple.pet_nondune.all.png",height=4000,width=4000)
as_tibble(data_pet[chr2 != chr1]) %>% 
  mutate(chr1 = gsub("Ha412HOChr","",chr1),
         chr2 = gsub("Ha412HOChr","",chr2)) %>%
  filter(chr1 > chr2) %>%
  mutate(nondune_links = case_when(nondune_links > 0.5 ~ 0.5,
                               TRUE ~ nondune_links)) %>%
  ggplot(.,aes()) + geom_tile(aes(x=bin1/1000000,y=bin2/1000000,fill=nondune_links)) +
  scale_fill_distiller(palette = "Spectral",name="HiC links") +
  theme_bw() + ylab("Chr6 Mbp") + xlab("Chr15 Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  ggtitle("Pet_nondune") +
  facet_grid(chr2~chr1,scales="free")
dev.off()

png("HiC/HiC_adaptTrim_NGM_homer10mq_simple.pet_comparison.all.png",height=4000,width=4000)
as_tibble(data_pet[chr2 != chr1]) %>% 
  mutate(chr1 = gsub("Ha412HOChr","",chr1),
         chr2 = gsub("Ha412HOChr","",chr2)) %>%
  filter(chr1 > chr2) %>%
  mutate(link_comparison = case_when(link_comparison > 0.2 ~ 0.2,
                                     link_comparison < -0.2 ~ -0.2,
                                   TRUE ~ link_comparison)) %>%
  ggplot(.,aes()) + geom_tile(aes(x=bin1/1000000,y=bin2/1000000,fill=link_comparison)) +
  scale_fill_distiller(palette = "Spectral",name="HiC link\ndifference") +
  theme_bw() + ylab("Chr6 Mbp") + xlab("Chr15 Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  ggtitle("Pet_dune_comparison") +
  facet_grid(chr2~chr1,scales="free")
dev.off()

png("HiC/HiC_adaptTrim_NGM_homer10mq_simple.petdunevann.all.png",height=4000,width=4000)
as_tibble(data_pet[chr2 == chr1 & bin1 != bin2]) %>%  
  mutate(dune_v_ann = case_when(dune_v_ann > 1 ~ 1,
                                dune_v_ann < -1 ~ -1,
                                TRUE ~ dune_v_ann)) %>%
  ggplot(.,aes()) + geom_tile(aes(x=bin1/1000000,y=bin2/1000000,fill=dune_v_ann)) +
  scale_fill_distiller(palette = "Spectral",name="HiC links",limits=c(-1,1)) +
  theme_bw() + ylab("Mbp") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  ggtitle("Pet_dune vs ANN") +
  facet_wrap(~chr1,scales="free")
dev.off()
png("HiC/HiC_adaptTrim_NGM_homer10mq_simple.petnondunevann.all.png",height=4000,width=4000)
as_tibble(data_pet[chr2 == chr1 & bin1 != bin2]) %>%  
  mutate(nondune_v_ann = case_when(nondune_v_ann > 1 ~ 1,
                                nondune_v_ann < -1 ~ -1,
                                TRUE ~ nondune_v_ann)) %>%
  ggplot(.,aes()) + geom_tile(aes(x=bin1/1000000,y=bin2/1000000,fill=nondune_v_ann)) +
  scale_fill_distiller(palette = "Spectral",name="HiC links",limits=c(-1,1)) +
  theme_bw() + ylab("Mbp") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  ggtitle("Pet_nondune vs ANN") +
  facet_wrap(~chr1,scales="free")
dev.off()
png("HiC/HiC_adaptTrim_NGM_homer10mq_simple.arginvvann.all.png",height=4000,width=4000)
as_tibble(data_arg[chr2 == chr1 & bin1 != bin2]) %>%  
  mutate(arginv_v_ann = case_when(arginv_v_ann > 1 ~ 1,
                                  arginv_v_ann < -1 ~ -1,
                                   TRUE ~ arginv_v_ann)) %>%
  ggplot(.,aes()) + geom_tile(aes(x=bin1/1000000,y=bin2/1000000,fill=arginv_v_ann)) +
  scale_fill_distiller(palette = "Spectral",name="HiC links",limits=c(-1,1)) +
  theme_bw() + ylab("Mbp") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  ggtitle("arg_inv vs ANN") +
  facet_wrap(~chr1,scales="free")
dev.off()
png("HiC/HiC_adaptTrim_NGM_homer10mq_simple.argnoninvvann.all.png",height=4000,width=4000)
as_tibble(data_arg[chr2 == chr1 & bin1 != bin2]) %>%  
  mutate(argnoninv_v_ann = case_when(argnoninv_v_ann > 1 ~ 1,
                                     argnoninv_v_ann < -1 ~ -1,
                                   TRUE ~ argnoninv_v_ann)) %>%
  ggplot(.,aes()) + geom_tile(aes(x=bin1/1000000,y=bin2/1000000,fill=argnoninv_v_ann)) +
  scale_fill_distiller(palette = "Spectral",name="HiC links",limits=c(-1,1)) +
  theme_bw() + ylab("Mbp") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  ggtitle("arg_noninv vs ANN") +
  facet_wrap(~chr1,scales="free")
dev.off()
