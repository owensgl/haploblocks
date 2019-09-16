library(tidyverse)
library(zoo)
library(SNPRelate)
library(grid)
library(gridExtra)
library(plotly)
library(scatterpie)
library(RColorBrewer)
#Visualize candidate inversions from lostruct

inversions <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt") %>% filter(spe != "niveus")
data_directory <- "/media/owens/Copper/wild_gwas"
window_size <- "100"
chr_prefix <- "Ha412HOChr"
het_script <- "/home/owens/bin/pop_gen/vcf2het.pl"
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T) 
pop_loc <- read_tsv("pop_loc_allnum.txt")
pop_loc %>% rename(population = pop) %>% inner_join(.,labels) %>% rename(sample = name) -> labels
info <- read_tsv("sample_info_apr_2018.tsv") %>% rename(sample = name)


read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt") %>%
  filter(spe != "niveus") %>%
  mutate(id = paste0(chr,".",mds)) -> inversion_locations
read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt") %>%
  filter(spe != "niveus") %>%
  mutate(id = paste0("Ha412HOChr",sprintf('%02d',chr),".",direction,mds)) %>%
  select(sv_name,id, spe) %>%
  inner_join(inversion_locations)-> inversion_locations

ld_plot_list = list()
pca_plot_list = list()
mds_plot_list = list()
het_plot_list = list()

#Top of page is n=1 and n=20
for (n in 1:nrow(inversions)){
  
  chosen_species_lower <- pull(inversions[n,1])
  chosen_species_upper <- pull(inversions[n,2])
  chosen_tag <- pull(inversions[n,3])
  chosen_chr_n <- sprintf("%02d",pull(inversions[n,4]))
  chosen_threshold <- pull(inversions[n,7])
  chosen_mds <- pull(inversions[n,5])
  chosen_direction <- pull(inversions[n,6])
  chosen_sv_name <- inversions$sv_name[n]
  
  chosen_start <- inversion_locations %>%
    filter(sv_name == chosen_sv_name) %>% pull(start) %>% min()
  chosen_end <- inversion_locations %>%
    filter(sv_name == chosen_sv_name) %>% pull(end) %>% max() 
  
  
  outline_color <- "black"
  #Decide on outline color
  if (chosen_species_lower == "annuus"){
    outline_color <- "#FFB31C"
  }else if (chosen_species_lower == "argophyllus"){
    outline_color <- "#3DA047"
  }else {
    outline_color <- "#447499"
  }
  #Make LD plot.Choose petiolaris subspecies that has more polymorphism
  if (chosen_species_lower == "petiolaris"){
    pet_species <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.petinvfreq.txt") %>%
      filter(chr == paste0("Ha412HOChr",chosen_chr_n),mds == paste0(chosen_direction,chosen_mds)) %>%
      mutate(maf = abs(freq - 0.5)) %>% arrange(maf) %>% head(1) %>% pull(species)
    if (pet_species == "PetFal"){
      chosen_label = "F"
    }else{
      chosen_label = "P"
    }
    between_ld <- read_tsv(paste0("/media/owens/Copper/wild_gwas/",tolower(pet_species),"/",
                                  chosen_species_upper, ".tranche90.snp.",tolower(pet_species),".90.bi.remappedHa412HO.thin100.maf5.Chr",
                                  chosen_chr_n,".windows.ld.gz")) %>%
      filter(win1 != win2)
    within_ld <- read_tsv(paste0("/media/owens/Copper/wild_gwas/",tolower(pet_species),"/",chosen_species_lower,".Ha412HOChr",chosen_chr_n,".",
                                 chosen_direction, chosen_mds,".withinhaplo.ld.txt.gz")) %>%
      rename(win1 = win2, win2 = win1) %>%
      filter(win1 != win2)
    

    ld_plot<-  rbind(between_ld,within_ld) %>%
      filter(n > 50) %>%
      ggplot(.,aes()) + 
      geom_tile(aes(x=win1/1000000,y=win2/1000000,fill=max_2_r2)) +
      annotate("segment",y=-4,x=floor(chosen_start/500000)/2,
               yend=-4,xend=floor(chosen_end/500000)/2,
               color="#9B2374",size=1) +
      annotate("segment",y=floor(chosen_start/500000)/2,x=-4,
               yend=floor(chosen_end/500000)/2,xend=-4,
               color="#9B2374",size=1) +
      scale_fill_viridis_c(limits=c(0,1),name="LD") +
      theme_linedraw() +ylab("Mbp") + xlab("Mbp") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            plot.margin = margin(0, 0, 2, 2, "pt")) +
      scale_x_continuous(limits = c(-4,max(between_ld$win2)/1000000), expand = c(0, 0)) +
      scale_y_continuous(limits = c(-4,max(between_ld$win2)/1000000), expand = c(0, 0)) +
      coord_cartesian(clip = 'off') +
      theme(plot.tag.position = c(0.15, 0.85),
            plot.tag=element_text(size=6,color="red")) 

    ld_plot_list[[n]] <- ld_plot
  }else{
    chosen_label = ""
    #For argophyllus and annuus
    between_ld <- read_tsv(paste0("/media/owens/Copper/wild_gwas/",chosen_species_lower,"/",
                                  chosen_species_upper, ".tranche90.snp.",chosen_tag,".90.bi.remappedHa412HO.thin100.maf5.Chr",
                                  chosen_chr_n,".windows.ld.gz")) %>%
      filter(win1 != win2)
    within_ld <- read_tsv(paste0("/media/owens/Copper/wild_gwas/",chosen_species_lower,"/",chosen_species_lower,".Ha412HOChr",chosen_chr_n,".",
                                 chosen_direction, chosen_mds,".withinhaplo.ld.txt.gz")) %>%
      rename(win1 = win2, win2 = win1) %>%
      filter(win1 != win2)
    ld_plot <- rbind(between_ld,within_ld) %>%
      filter(n > 50) %>%
      ggplot(.,aes()) + 
      geom_tile(aes(x=win1/1000000,y=win2/1000000,fill=max_2_r2)) +
      annotate("segment",y=-4,x=floor(chosen_start/500000)/2,
               yend=-4,xend=floor(chosen_end/500000)/2,
               color="#9B2374",size=1) +
      annotate("segment",y=floor(chosen_start/500000)/2,x=-4,
               yend=floor(chosen_end/500000)/2,xend=-4,
               color="#9B2374",size=1) +
      scale_fill_viridis_c(limits=c(0,1),name="LD") +
      theme_linedraw() +ylab("Mbp") + xlab("Mbp") +

      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            plot.tag.position = c(.05, 0.12),
            plot.margin = margin(0, 0, 2, 2, "pt")) +
      scale_x_continuous(limits = c(-4,max(between_ld$win2)/1000000), expand = c(0, 0)) +
      scale_y_continuous(limits = c(-4,max(between_ld$win2)/1000000), expand = c(0, 0)) +
      coord_cartesian(clip = 'off') +
      theme(plot.tag.position = c(0.15, 0.85),
            plot.tag=element_text(size=6,color="red")) 
    
    ld_plot_list[[n]] <- ld_plot
  }
    
    
  #For petiolaris, we have to pick a single group (PetPF, PetPet, or PetFal) for each of the synchronized regions
  if (chosen_species_lower == "petiolaris"){
    guide <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.synchronizePet.txt")
    old_inversions <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt")
    guide %>%
      filter(Chr == paste0("Ha412HOChr",chosen_chr_n), synchronized==paste0(chosen_direction,chosen_mds)) -> filtered_guide
    if(!is.na(filtered_guide$PetPF[1])){
      #It's present in PetPF
      old_inversions %>% filter(spe == "petiolaris",chr== as.numeric(chosen_chr_n)) %>%
        mutate(mds_together = paste0(direction,mds)) %>%
        filter(mds_together== filtered_guide$PetPF) -> new_info
      chosen_species_lower <- pull(new_info[1,1])
      chosen_species_upper <- pull(new_info[1,2])
      chosen_tag <- pull(new_info[1,3])
      chosen_chr_n <- sprintf("%02d",pull(new_info[1,4]))
      chosen_threshold <- pull(new_info[1,7])
      chosen_mds <- pull(new_info[1,5])
      chosen_direction <- pull(new_info[1,6])
      chosen_label <- "PF"
    }else if (!is.na(filtered_guide$PetPet[1])){
      old_inversions %>% filter(spe == "petpet",chr== as.numeric(chosen_chr_n)) %>%
        mutate(mds_together = paste0(direction,mds)) %>%
        filter(mds_together== filtered_guide$PetPet) -> new_info
      chosen_species_lower <- pull(new_info[1,1])
      chosen_species_upper <- pull(new_info[1,2])
      chosen_tag <- pull(new_info[1,3])
      chosen_chr_n <- sprintf("%02d",pull(new_info[1,4]))
      chosen_threshold <- pull(new_info[1,7])
      chosen_mds <- pull(new_info[1,5])
      chosen_direction <- pull(new_info[1,6])
      chosen_label <- "P"
      
    }else{
      old_inversions %>% filter(spe == "petfal",chr== as.numeric(chosen_chr_n)) %>%
        mutate(mds_together = paste0(direction,mds)) %>%
        filter(mds_together== filtered_guide$PetFal) -> new_info
      chosen_species_lower <- pull(new_info[1,1])
      chosen_species_upper <- pull(new_info[1,2])
      chosen_tag <- pull(new_info[1,3])
      chosen_chr_n <- sprintf("%02d",pull(new_info[1,4]))
      chosen_threshold <- pull(new_info[1,7])
      chosen_mds <- pull(new_info[1,5])
      chosen_direction <- pull(new_info[1,6])
      chosen_label <- "F"
      
    }
    
  }
  win.regions <- readRDS(paste(
    "MDS_plots/Ha412HO/",
    chosen_species_lower,
    "/",
    chosen_species_upper,
    ".tranche90.snp.",
    chosen_tag,
    ".90.bi.remappedHa412HO.Chr",
    chosen_chr_n,
    ".",
    window_size,
    ".windows.rds",sep=""))
  
  win.regions %>%
    select(chrom,start,end,n,mid,paste0("mds",sprintf("%02d",chosen_mds))) %>%
    rename(mds_chosen = paste0("mds",sprintf("%02d",chosen_mds)))-> win.regions
  
  if (chosen_direction == "neg"){
    win.regions <- win.regions %>%
      mutate(outlier = case_when(mds_chosen < chosen_threshold ~ "1",
                                 TRUE ~ "0")) 
  }else{
    win.regions <- win.regions %>%
      mutate(outlier = case_when(mds_chosen > chosen_threshold ~ "1",
                                 TRUE ~ "0")) 
  }
  
  #Select windows surrounded by other outliers.
  max_distance_surround <- 20
  win.regions <- win.regions %>% mutate(lead_max = rollmax(x = outlier, max_distance_surround, align = "left", fill = NA),
                                        lag_max = rollmax(x = outlier, max_distance_surround, align = "right", fill = NA)) %>%
    mutate(filled_outlier = case_when(outlier == 1 ~ 1,
                                      lead_max == 1 & lag_max == 1 ~ 1,
                                      TRUE ~ 0)) 
  #Plot mds scores
  mds_plot <- win.regions %>%
    ggplot(.,aes(x=mid/1000000,y=mds_chosen,color=as.factor(filled_outlier))) + 
    annotate(geom="segment",x=0,xend=max(win.regions$mid)/1000000,y=0,yend=0,color="grey") +
    geom_point(size=0.1) +
    theme_bw() + scale_color_manual(values=c("black",outline_color)) +
    ylab(paste("mds",chosen_chr_n,sep="")) +
    xlab("Mb") +
    theme(legend.position = "none") +
    theme(panel.grid.major.y = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          #axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          text = element_text(size=8),
          panel.border = element_rect(colour = "white", fill=NA, size=0),
          plot.margin = margin(0, 10, 0, 0, "pt")) +
    annotate(geom="segment",x=0,xend=0,y=max(win.regions$mds_chosen)+0.1,yend=min(win.regions$mds_chosen)-0.1) +
    annotate(geom="text",x=-10,,y=((max(win.regions$mds_chosen)-min(win.regions$mds_chosen))/2)+ min(win.regions$mds_chosen),label=chosen_sv_name,
             angle = 90,size=1.75) +
    theme(plot.tag.position = c(0.03, 0.92),
          plot.tag=element_text(size=8,face="bold")) 
  
  if (n == 1){
    mds_plot <- mds_plot + labs(tag="a")
  }
  mds_plot_list[[n]] <- mds_plot
  
  win.regions$ungapped_n <- 1:nrow(win.regions)
  win.regions <- win.regions %>% 
    filter(filled_outlier == 1) %>% mutate(gap = ungapped_n - lag(ungapped_n))
  
  win.regions$group <- "1"
  current_group = 1
  for (i in 2:nrow(win.regions)){
    if (win.regions$gap[i] > 1){
      current_group = current_group + 1
    }
    win.regions$group[i] <- current_group
  }
  win.regions %>% 
    group_by(group) %>%
    summarize(start = min(start),end=max(end)) %>%
    mutate(chr=paste(chr_prefix,chosen_chr_n,sep="")) -> region_selected
  
  #Now return to the full species set specifically for petiolaris (which is plotting the MDS from a subsample)
  chosen_species_lower <- pull(inversions[n,1])
  chosen_species_upper <- pull(inversions[n,2])
  chosen_tag <- pull(inversions[n,3])
  chosen_chr_n <- sprintf("%02d",pull(inversions[n,4]))
  chosen_threshold <- pull(inversions[n,7])
  chosen_mds <- pull(inversions[n,5])
  chosen_direction <- pull(inversions[n,6])
  
  if (chosen_species_lower == "petiolaris"){
    chosen_tag <- "petPF"
  }
  
  #Check if the PCA is already saved.
  if(!file.exists(paste0(
    "MDS_outliers/Ha412HO/",
    chosen_species_lower,
    "/",
    "Ha412HO_inv.v3.pcasites.Ha412HOChr",
    chosen_chr_n,
    ".",
    chosen_direction,
    chosen_mds,
    ".pca.txt"))){
    
    vcf_in <- paste(
      data_directory,
      "/",
      chosen_species_lower,
      "/",
      chosen_species_upper,
      ".tranche90.snp.",
      chosen_tag,
      ".90.bi.remappedHa412HO.vcf.gz",
      sep="")
    vcf_out <- paste(
      data_directory,
      "/",
      chosen_species_lower,
      "/",
      chosen_species_upper,
      ".tranche90.snp.",
      chosen_tag,
      ".90.bi.remappedHa412HO.tmp.vcf.gz",
      sep="")
    gds_out <- paste(
      data_directory,
      "/",
      chosen_species_lower,
      "/",
      chosen_species_upper,
      ".tranche90.snp.",
      chosen_tag,
      ".90.bi.remappedHa412HO.tmp.gds",
      sep="")
    region_string <- paste(region_selected$chr[1],":",region_selected$start[1],"-",region_selected$end[1],sep="")
    if (nrow(region_selected)> 1){
      for (i in 2:nrow(region_selected)){
        region_string <- paste(region_string,",",region_selected$chr[i],":",region_selected$start[i],"-",region_selected$end[i],sep="")
      }
    }
    
    
    system(paste("bcftools view -O z -r",region_string,vcf_in,">",vcf_out))
    system(paste("tabix -p vcf",vcf_out))
    
    system(paste("zcat ",vcf_out," | perl ",het_script," > tmp.het.txt",sep=""))
    vcf_het <- read_tsv("tmp.het.txt")
    snpgdsVCF2GDS(vcf_out,gds_out, method="biallelic.only",ignore.chr.prefix="Ha412HOChr")
    genofile <- snpgdsOpen(gds_out)
    pca <- snpgdsPCA(genofile, num.thread=2)
    tab <- data.frame(sample = pca$sample.id,
                      EV1 = pca$eigenvect[,1],    # the first eigenvector
                      EV2 = pca$eigenvect[,2],    # the second eigenvector
                      stringsAsFactors = FALSE)
    snpgdsClose(genofile)
    samplelist <- as.vector(tab$sample)
    inv_genotypes <- read_tsv(paste(
      "MDS_outliers/Ha412HO/",
      chosen_species_lower,
      "/",
      "Ha412HO_inv.v3.pcasites.Ha412HOChr",
      chosen_chr_n,
      ".",
      chosen_direction,
      chosen_mds,
      ".genotypes.txt",sep=""))
    
    
    ##Check to see if we have to flip the first PCA, since it randomly assigns it each time it is run and needs to be consistent with the genotypes file
    tab %>% 
      inner_join(inv_genotypes) %>%
      group_by(triangle_genotype) %>%
      summarize(mean_EV1 = mean(EV1)) -> EV1_summary
    if (EV1_summary$mean_EV1[1] > EV1_summary$mean_EV1[3]){
      tab$EV1 <- -tab$EV1
    }
    tab %>% inner_join(.,vcf_het) -> tab
    ##Save the PCA so it doesn't have to be rerun
    write_tsv(tab,paste0(
      "MDS_outliers/Ha412HO/",
      chosen_species_lower,
      "/",
      "Ha412HO_inv.v3.pcasites.Ha412HOChr",
      chosen_chr_n,
      ".",
      chosen_direction,
      chosen_mds,
      ".pca.txt"))
  }else{
    tab <- read_tsv(paste0(
      "MDS_outliers/Ha412HO/",
      chosen_species_lower,
      "/",
      "Ha412HO_inv.v3.pcasites.Ha412HOChr",
      chosen_chr_n,
      ".",
      chosen_direction,
      chosen_mds,
      ".pca.txt"))
  }
  inv_genotypes <- read_tsv(paste(
    "MDS_outliers/Ha412HO/",
    chosen_species_lower,
    "/",
    "Ha412HO_inv.v3.pcasites.Ha412HOChr",
    chosen_chr_n,
    ".",
    chosen_direction,
    chosen_mds,
    ".genotypes.txt",sep=""))
  
  pca_plot <- tab %>% 
    inner_join(inv_genotypes) %>%
    ggplot(.,aes(x=EV1,y=EV2)) + 
    geom_point(aes(color=as.factor(triangle_genotype)),size=0.15) +
    theme_bw() + 
    scale_color_manual(values=c("black","dark grey","light grey")) +
    theme(legend.position="none") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x= element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.title.y= element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          panel.border = element_rect(colour = outline_color, fill=NA, size=0.5),
          plot.margin = margin(2, 0, 2, 0, "pt"),
          text = element_text(size=6)) +
    ylab("PC2") + xlab("PC1") +
    theme(plot.tag.position = c(0.01, 0.94))
  
  

  pca_plot_list[[n]] <- pca_plot
  
  
  het_plot <- tab %>% 
    inner_join(inv_genotypes) %>%
    mutate(new_triangle_genotype = case_when(triangle_genotype == 0 ~ '0/0',
                                              triangle_genotype == 1 ~ '0/1',
                                              triangle_genotype == 2 ~ '1/1')) %>%
    ggplot(.,aes(x=new_triangle_genotype,y=percent_het)) + 
    geom_boxplot(aes(fill=new_triangle_genotype),size=0.1,outlier.size=0.2) +
    ylab(expression(H[o])) +
    theme_bw() + 
    scale_fill_manual(values=c("black","dark grey","light grey")) +
    theme(legend.position="none") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank(),
          text = element_text(size=6),
          panel.border = element_rect(colour = outline_color, fill=NA, size=0.5),
          plot.margin = margin(1, 1, 1, 5, "pt")) +
    theme(plot.tag.position = c(0.01, 0.94))
  

  het_plot_list[[n]] <- het_plot
  
  
  
 
}
#Make wave example plot
win.regions <- readRDS(paste(
  "MDS_plots/Ha412HO/",
  "annuus",
  "/",
  "Annuus",
  ".tranche90.snp.",
  "env",
  ".90.bi.remappedHa412HO.Chr",
  "07",
  ".",
  "100",
  ".windows.rds",sep=""))

wave_plot <- win.regions %>%
  ggplot(.,aes(x=mid/1000000,y=mds03)) + 
  annotate(geom="segment",x=0,xend=max(win.regions$mid)/1000000,y=0,yend=0,color="grey") +
  geom_point(size=0.1) +
  theme_bw() + scale_color_manual(values=c("black",outline_color)) +
  ylab(paste("mds",chosen_chr_n,sep="")) +
  xlab("Mb") +
  theme(legend.position = "none") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        #axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        text = element_text(size=8),
        panel.border = element_rect(colour = "white", fill=NA, size=0),
        plot.margin = margin(0, 10, 0, 10, "pt"),
        plot.tag=element_text(size=8,face="bold")) +
  annotate(geom="segment",x=0,xend=0,y=max(win.regions$mds03)+0.1,yend=min(win.regions$mds03)-0.1) +
  annotate(geom="text",x=-10,,y=((max(win.regions$mds03)-min(win.regions$mds03))/2)+ min(win.regions$mds03),label="mds03",
           angle = 90,size=1.75) +
  labs(tag="b")


# 
# 
# pdf("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.figureplots.ann.pdf",height=9.7,width=7.2)
# 
# grid.arrange(mds_plot_list[[1]], pca_plot_list[[1]], het_plot_list[[1]],
#              mds_plot_list[[2]], pca_plot_list[[2]], het_plot_list[[2]],
#              mds_plot_list[[3]], pca_plot_list[[3]], het_plot_list[[3]],
#              mds_plot_list[[4]], pca_plot_list[[4]], het_plot_list[[4]],
#              mds_plot_list[[5]], pca_plot_list[[5]], het_plot_list[[5]],
#              mds_plot_list[[6]], pca_plot_list[[6]], het_plot_list[[6]],
#              mds_plot_list[[7]], pca_plot_list[[7]], het_plot_list[[7]],
#              mds_plot_list[[8]], pca_plot_list[[8]], het_plot_list[[8]],
#              mds_plot_list[[9]], pca_plot_list[[9]], het_plot_list[[9]],
#              mds_plot_list[[10]], pca_plot_list[[10]], het_plot_list[[10]],
#              mds_plot_list[[11]], pca_plot_list[[11]], het_plot_list[[11]],
#              layout_matrix = rbind(c(1,1,1,2,3,3),
#                                    c(4,4,4,5,6,6),
#                                    c(7,7,7,8,9,9),
#                                    c(10,10,10,11,12,12),
#                                    c(13,13,13,14,15,15),
#                                    c(16,16,16,17,18,18),
#                                    c(19,19,19,20,21,21),
#                                    c(22,22,22,23,24,24),
#                                    c(25,25,25,26,27,27),
#                                    c(28,28,28,29,30,30),
#                                    c(31,31,31,32,33,33)))
# dev.off()
# 
# pdf("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.figureplots.annarg.pdf",height=7,width=7.2)
# 
# grid.arrange(mds_plot_list[[1]], pca_plot_list[[1]], het_plot_list[[1]],
#              mds_plot_list[[2]], pca_plot_list[[2]], het_plot_list[[2]],
#              mds_plot_list[[3]], pca_plot_list[[3]], het_plot_list[[3]],
#              mds_plot_list[[4]], pca_plot_list[[4]], het_plot_list[[4]],
#              mds_plot_list[[5]], pca_plot_list[[5]], het_plot_list[[5]],
#              mds_plot_list[[6]], pca_plot_list[[6]], het_plot_list[[6]],
#              mds_plot_list[[7]], pca_plot_list[[7]], het_plot_list[[7]],
#              mds_plot_list[[8]], pca_plot_list[[8]], het_plot_list[[8]],
#              mds_plot_list[[9]], pca_plot_list[[9]], het_plot_list[[9]],
#              mds_plot_list[[10]], pca_plot_list[[10]], het_plot_list[[10]],
#              mds_plot_list[[11]], pca_plot_list[[11]], het_plot_list[[11]],
#              mds_plot_list[[12]], pca_plot_list[[12]], het_plot_list[[12]],
#              mds_plot_list[[13]], pca_plot_list[[13]], het_plot_list[[13]],
#              mds_plot_list[[14]], pca_plot_list[[14]], het_plot_list[[14]],
#              mds_plot_list[[15]], pca_plot_list[[15]], het_plot_list[[15]],
#              mds_plot_list[[16]], pca_plot_list[[16]], het_plot_list[[16]],
#              mds_plot_list[[17]], pca_plot_list[[17]], het_plot_list[[17]],
#              layout_matrix = rbind(c(1,1,1,2,3,3,28,28,28,29,30,30),
#                                    c(4,4,4,5,6,6,31,31,31,32,33,33),
#                                    c(7,7,7,8,9,9,34,34,34,35,36,36),
#                                    c(10,10,10,11,12,12,37,37,37,38,39,39),
#                                    c(13,13,13,14,15,15,40,40,40,41,42,42),
#                                    c(16,16,16,17,18,18,43,43,43,44,45,45),
#                                    c(19,19,19,20,21,21,46,46,46,47,48,48),
#                                    c(22,22,22,23,24,24,49,49,49,50,51,51),
#                                    c(25,25,25,26,27,27,NA,NA,NA,NA,NA,NA)))
# 
# dev.off()
# 
# 
# pdf("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.figureplots.arg.pdf",height=6,width=5)
# 
# grid.arrange(mds_plot_list[[12]], pca_plot_list[[12]], het_plot_list[[12]],
#              mds_plot_list[[13]], pca_plot_list[[13]], het_plot_list[[13]],
#              mds_plot_list[[14]], pca_plot_list[[14]], het_plot_list[[14]],
#              mds_plot_list[[15]], pca_plot_list[[15]], het_plot_list[[15]],
#              mds_plot_list[[16]], pca_plot_list[[16]], het_plot_list[[16]],
#              mds_plot_list[[17]], pca_plot_list[[17]], het_plot_list[[17]],
#              layout_matrix = rbind(c(1,1,1,2,3,3),
#                                    c(4,4,4,5,6,6),
#                                    c(7,7,7,8,9,9),
#                                    c(10,10,10,11,12,12),
#                                    c(13,13,13,14,15,15),
#                                    c(16,16,16,17,18,18)))
# 
# 
# dev.off()
# 
# pdf("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.figureplots.pet.pdf",height=8.75,width=7.2)
# 
# grid.arrange(mds_plot_list[[18]], pca_plot_list[[18]], het_plot_list[[18]],
#              mds_plot_list[[19]], pca_plot_list[[19]], het_plot_list[[19]],
#              mds_plot_list[[20]], pca_plot_list[[20]], het_plot_list[[20]],
#              mds_plot_list[[21]], pca_plot_list[[21]], het_plot_list[[21]],
#              mds_plot_list[[22]], pca_plot_list[[22]], het_plot_list[[22]],
#              mds_plot_list[[23]], pca_plot_list[[23]], het_plot_list[[23]],
#              mds_plot_list[[24]], pca_plot_list[[24]], het_plot_list[[24]],
#              mds_plot_list[[25]], pca_plot_list[[25]], het_plot_list[[25]],
#              mds_plot_list[[26]], pca_plot_list[[26]], het_plot_list[[26]],
#              mds_plot_list[[27]], pca_plot_list[[27]], het_plot_list[[27]],
#              mds_plot_list[[28]], pca_plot_list[[28]], het_plot_list[[28]],
#              mds_plot_list[[29]], pca_plot_list[[29]], het_plot_list[[29]],
#              mds_plot_list[[30]], pca_plot_list[[30]], het_plot_list[[30]],
#              mds_plot_list[[31]], pca_plot_list[[31]], het_plot_list[[31]],
#              mds_plot_list[[32]], pca_plot_list[[32]], het_plot_list[[32]],
#              mds_plot_list[[33]], pca_plot_list[[33]], het_plot_list[[33]],
#              mds_plot_list[[34]], pca_plot_list[[34]], het_plot_list[[34]],
#              mds_plot_list[[35]], pca_plot_list[[35]], het_plot_list[[35]],
#              mds_plot_list[[36]], pca_plot_list[[36]], het_plot_list[[36]],
#              mds_plot_list[[37]], pca_plot_list[[37]], het_plot_list[[37]],
#              layout_matrix = rbind(c(1,1,1,2,3,3,31,31,31,32,33,33),
#                                    c(4,4,4,5,6,6,34,34,34,35,36,36),
#                                    c(7,7,7,8,9,9,37,37,37,38,39,39),
#                                    c(10,10,10,11,12,12,40,40,40,41,42,42),
#                                    c(13,13,13,14,15,15,43,43,43,44,45,45),
#                                    c(16,16,16,17,18,18,46,46,46,47,48,48),
#                                    c(19,19,19,20,21,21,49,49,49,50,51,51),
#                                    c(22,22,22,23,24,24,52,52,52,53,54,54),
#                                    c(25,25,25,26,27,27,55,55,55,56,57,57),
#                                    c(28,28,28,29,30,30,58,58,58,59,60,60)))
# dev.off()

pdf("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.figureplots.all.1page.pdf",height=9.7,width=7.2)

#Make layout grid
layout <- matrix(, nrow = 0, ncol = 16)
for (x in 0:18){
  base <- x*4
  extra <- (x+19)*4
  if( x < 18){
    row <- c(base, base, base, base+1, base+2,base+2,base+3, base+3,extra,extra,extra,extra+1,extra+2,extra+2,extra+3,extra+3)
    
  }else{
    row <- c(base, base, base, base+1, base+2, base+2,base+3,base+3,148,148,148,148,148,148,148,148)
    
  }
  layout <- rbind(layout,row)
}

grid.arrange(mds_plot_list[[1]], pca_plot_list[[1]], het_plot_list[[1]], ld_plot_list[[1]],
             mds_plot_list[[2]], pca_plot_list[[2]], het_plot_list[[2]], ld_plot_list[[2]],
             mds_plot_list[[3]], pca_plot_list[[3]], het_plot_list[[3]], ld_plot_list[[3]],
             mds_plot_list[[4]], pca_plot_list[[4]], het_plot_list[[4]], ld_plot_list[[4]],
             mds_plot_list[[5]], pca_plot_list[[5]], het_plot_list[[5]], ld_plot_list[[5]],
             mds_plot_list[[6]], pca_plot_list[[6]], het_plot_list[[6]], ld_plot_list[[6]],
             mds_plot_list[[7]], pca_plot_list[[7]], het_plot_list[[7]], ld_plot_list[[7]],
             mds_plot_list[[8]], pca_plot_list[[8]], het_plot_list[[8]], ld_plot_list[[8]],
             mds_plot_list[[9]], pca_plot_list[[9]], het_plot_list[[9]], ld_plot_list[[9]],
             mds_plot_list[[10]], pca_plot_list[[10]], het_plot_list[[10]], ld_plot_list[[10]],
             mds_plot_list[[11]], pca_plot_list[[11]], het_plot_list[[11]], ld_plot_list[[11]],
             mds_plot_list[[12]], pca_plot_list[[12]], het_plot_list[[12]], ld_plot_list[[12]],
             mds_plot_list[[13]], pca_plot_list[[13]], het_plot_list[[13]], ld_plot_list[[13]],
             mds_plot_list[[14]], pca_plot_list[[14]], het_plot_list[[14]], ld_plot_list[[14]],
             mds_plot_list[[15]], pca_plot_list[[15]], het_plot_list[[15]], ld_plot_list[[15]],
             mds_plot_list[[16]], pca_plot_list[[16]], het_plot_list[[16]], ld_plot_list[[16]],
             mds_plot_list[[17]], pca_plot_list[[17]], het_plot_list[[17]], ld_plot_list[[17]],
             mds_plot_list[[18]], pca_plot_list[[18]], het_plot_list[[18]], ld_plot_list[[18]],
             mds_plot_list[[19]], pca_plot_list[[19]], het_plot_list[[19]], ld_plot_list[[19]],
             mds_plot_list[[20]], pca_plot_list[[20]], het_plot_list[[20]], ld_plot_list[[20]],
             mds_plot_list[[21]], pca_plot_list[[21]], het_plot_list[[21]], ld_plot_list[[21]],
             mds_plot_list[[22]], pca_plot_list[[22]], het_plot_list[[22]], ld_plot_list[[22]],
             mds_plot_list[[23]], pca_plot_list[[23]], het_plot_list[[23]], ld_plot_list[[23]],
             mds_plot_list[[24]], pca_plot_list[[24]], het_plot_list[[24]], ld_plot_list[[24]],
             mds_plot_list[[25]], pca_plot_list[[25]], het_plot_list[[25]], ld_plot_list[[25]],
             mds_plot_list[[26]], pca_plot_list[[26]], het_plot_list[[26]], ld_plot_list[[26]],
             mds_plot_list[[27]], pca_plot_list[[27]], het_plot_list[[27]], ld_plot_list[[27]],
             mds_plot_list[[28]], pca_plot_list[[28]], het_plot_list[[28]], ld_plot_list[[28]],
             mds_plot_list[[29]], pca_plot_list[[29]], het_plot_list[[29]], ld_plot_list[[29]],
             mds_plot_list[[30]], pca_plot_list[[30]], het_plot_list[[30]], ld_plot_list[[30]],
             mds_plot_list[[31]], pca_plot_list[[31]], het_plot_list[[31]], ld_plot_list[[31]],
             mds_plot_list[[32]], pca_plot_list[[32]], het_plot_list[[32]], ld_plot_list[[32]],
             mds_plot_list[[33]], pca_plot_list[[33]], het_plot_list[[33]], ld_plot_list[[33]],
             mds_plot_list[[34]], pca_plot_list[[34]], het_plot_list[[34]], ld_plot_list[[34]],
             mds_plot_list[[35]], pca_plot_list[[35]], het_plot_list[[35]], ld_plot_list[[35]],
             mds_plot_list[[36]], pca_plot_list[[36]], het_plot_list[[36]], ld_plot_list[[36]],
             mds_plot_list[[37]], pca_plot_list[[37]], het_plot_list[[37]], ld_plot_list[[37]],
             wave_plot,
             layout_matrix = layout)
dev.off()



layout_matrix = rbind(c(1,1,1,2,3,3,58,58,58,59,60,60),
                      c(4,4,4,5,6,6,61,61,61,62,63,63),
                      c(7,7,7,8,9,9,64,64,64,65,66,66),
                      c(10,10,10,11,12,12,67,67,67,68,69,69),
                      c(13,13,13,14,15,15,70,70,70,71,72,72),
                      c(16,16,16,17,18,18,73,73,73,74,75,75),
                      c(19,19,19,20,21,21,76,76,76,77,78,78),
                      c(22,22,22,23,24,24,79,79,79,80,81,81),
                      c(25,25,25,26,27,27,82,82,82,83,84,84),
                      c(28,28,28,29,30,30,85,85,85,86,87,87),
                      c(31,31,31,32,33,33,88,88,88,89,90,90),
                      c(34,34,34,35,36,36,91,91,91,92,93,93),
                      c(37,37,37,38,39,39,94,94,94,95,96,96),
                      c(40,40,40,41,42,42,97,97,97,98,99,99),
                      c(43,43,43,44,45,45,100,100,100,101,102,102),
                      c(46,46,46,47,48,48,103,103,103,104,105,105),
                      c(49,49,49,50,51,51,106,106,106,107,108,108),
                      c(52,52,52,53,54,54,109,109,109,110,111,111),
                      c(55,55,55,56,57,57,NA,NA,NA,NA,NA,NA)))

             