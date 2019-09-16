#Comparing the GEA with SV in correlation matrix and without them
library(tidyverse)
library(ggrepel)



gff <- read_tsv("/home/owens/ref/Han412-HO_gene_filt.gff3.gz",comment="#",
                col_names = c("chr","spacer","type","start","end","spacer2","spacer3","spacer4","info")) %>%
  filter(type == "gene") %>%
  separate(info,c("id","name"),";") %>%
  mutate(id = gsub("ID=gene:","",id))%>%
  select(chr,start,end,id)

space_around_genes = 1000;


species_list <- c("annuus","argophyllus","petfal","petpet")
capital_species_list <- c("Annuus","Argophyllus","Petiolaris","Petiolaris")
chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")
min_freq <- 0.03
for (x in 1:4){
  chosen_species <- species_list[x]
  chosen_capital_species <- capital_species_list[x]
  climate_column_names <- colnames(read_tsv(paste0("gea/column.names.",chosen_species,".varout.txt")))
  soil_column_names <- colnames(read_tsv("/media/owens/Copper/wild_gwas/env_associations/soil/header_info.txt"))
  
  climate_associations_with_sv <- read_tsv(paste0("/media/owens/Copper/wild_gwas/env_associations/climate/snps_HA412/with_inv_covariates/",
                                      chosen_species,"/var_out_",chosen_species,"_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.txt.gz"),
                               col_names = climate_column_names) %>%
    select(-chrom,-pos) %>%
    separate(snp_id, c("chr","pos"),"__") %>%
    mutate(pos = as.numeric(pos)) 
  climate_associations_without_sv <- read_tsv(paste0("/media/owens/Copper/wild_gwas/env_associations/climate/snps_HA412/without_inv_covariates/",
                                              chosen_species,"/var_out_",chosen_species,"_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD_noinv.txt.gz"),
                                       col_names = climate_column_names) %>%
    select(-chrom,-pos) %>%
    separate(snp_id, c("chr","pos"),"__") %>%
    mutate(pos = as.numeric(pos)) 
  
  soil_associations_with_sv <- read_delim(paste0("/media/owens/Copper/wild_gwas/env_associations/soil/snps_HA412/with_inv_covariates/",
                                               chosen_species,"/",chosen_species,"_HA412_reg_matrix_all_covariates_baypass.BF.txt.gz"),
                                       col_names = soil_column_names,delim=" ") %>%
    select(-chr,-pos) %>%
    separate(chr_pos, c("chr","pos"),"__") %>%
    mutate(pos = as.numeric(pos)) 
  soil_associations_without_sv <- read_delim(paste0("/media/owens/Copper/wild_gwas/env_associations/soil/snps_HA412/without_inv_covariates/",
                                                 chosen_species,"/",chosen_species,"_HA412_invfree_matrix_all_covariates_baypass.BF.txt.gz"),
                                          col_names = soil_column_names,delim=" ") %>%
    select(-chr,-pos) %>%
    separate(chr_pos, c("chr","pos"),"__") %>%
    mutate(pos = as.numeric(pos)) 
  
  cum_location <- chr_lengths %>% 
    select(chr,end) %>%
    rename(chr_len = end) %>%
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(rbind(climate_associations_with_sv %>% select(chr, pos),
                    soil_associations_with_sv %>% select(chr, pos)) %>%
                unique(), by=c("chr"="chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(chr, pos) %>%
    mutate( poscum=pos+tot)
  
  inner_join(climate_associations_without_sv %>%
    gather(.,variable,bayes_without_sv,latitude_e:RH_e)  %>%
    filter(variable != "PAS_e"),
    climate_associations_with_sv %>%
    gather(.,variable,bayes_with_sv,latitude_e:RH_e)  %>%
    filter(variable != "PAS_e")) -> climate_associations_combined
  
  inner_join(soil_associations_without_sv %>%
               gather(.,variable,bayes_without_sv,OM:SOL_SALTS),
             soil_associations_with_sv %>%
               gather(.,variable,bayes_with_sv,OM:SOL_SALTS)) -> soil_associations_combined
  
  axisdf = cum_location %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )
  
  #Load up inversions and their scores.
  
  climate_associations <- read_tsv(paste0("/media/owens/Copper/wild_gwas/env_associations/climate/sv_HA412/",chosen_species,
                                      "/varout_",chosen_species,"_HA412_remapped_Baypass_noinv_matrix_25_vars.tab")) %>%
    mutate(full_id = inversion_ID) %>%
    mutate(inversion_ID = gsub("Ha412HOChr","",inversion_ID)) %>%
    mutate(inversion_ID = gsub("pos","",inversion_ID)) %>%
    mutate(inversion_ID = gsub("neg","",inversion_ID)) %>%
    mutate(inversion_ID = gsub("syn","",inversion_ID)) %>%
    gather(variable,bayesfactor,colnames(.)[2]:colnames(.)[ncol(.)-1])
  
  soil_associations <-read_tsv(paste0("/media/owens/Copper/wild_gwas/env_associations/soil/sv_HA412/baypass_inversions_soil_data_bf.txt")) %>%
    filter(species == chosen_species) %>%
      mutate(full_id = paste0(chr,"_",mds)) %>%
      mutate(inversion_ID = full_id) %>%
      mutate(inversion_ID = gsub("Ha412HOChr","",inversion_ID)) %>%
      mutate(inversion_ID = gsub("pos","",inversion_ID)) %>%
      mutate(inversion_ID = gsub("neg","",inversion_ID)) %>%
      mutate(inversion_ID = gsub("syn","",inversion_ID)) %>%
      select(-ID,-species,-mds) %>%
      gather(variable,bayesfactor,colnames(.)[2]:colnames(.)[ncol(.)-2]) %>%
    select(-chr)
  
  if (chosen_species == "petpet" | chosen_species == "petfal"){
    #Remove SVs below a minimum MAF in petiolaris
    inv_frequencies <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.petinvfreq.txt") %>%
      filter(tolower(species) == chosen_species) %>%
      mutate(chr = gsub("Ha412HOChr","",chr),
             pos=gsub("syn","",mds)) %>%
      select(-mds,-species) %>%
      mutate(inversion_ID = paste0(chr,"_",pos))
    inner_join(climate_associations, inv_frequencies) %>%
      filter(freq > min_freq,freq < (1-min_freq)) %>%
      select(-freq) -> climate_associations
    inner_join(soil_associations, inv_frequencies) %>%
      filter(freq > min_freq,freq < (1-min_freq)) %>%
      select(-freq) -> soil_associations
  }
  if (chosen_species == "annuus"){
    species_full_name <- "H. annuus"
  }else if (chosen_species == "argophyllus"){
    species_full_name <- "H. argophyllus"
  }else if (chosen_species == "petpet"){
    species_full_name <- "H. petiolaris petiolaris"
  }else if (chosen_species == "petfal"){
    species_full_name <- "H. petiolaris fallax"
  }
  #Load up inversion locations
  inv_locations <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt") %>%
    mutate(full_id = paste0(chr,"_",mds)) %>%
    filter(species == chosen_capital_species)
  
  chr_lengths %>% 
    select(chr,end) %>%
    rename(chr_len = end) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(inv_locations, ., by=c("chr"="chr")) %>%
    # Add a cumulative position of each SNP
    arrange(chr, start) %>%
    mutate( startcum=start+tot,endcum=end+tot) %>%
    separate(chr, c("ref","chr_n"), "Chr") %>%
    mutate(chr_n = as.numeric(chr_n)) -> inv_locations
  
 
  
  unique(climate_associations$variable) -> climate_variables
  unique(soil_associations$variable) -> soil_variables
  
  
  
  ###########
  #Getting genes with >20 dB
  ###########
  soil_variables <- unique(soil_associations_combined$variable)
  climate_variables <- unique(climate_associations_combined$variable)
  bayes_20_genes <- tibble(chr=character(),start=numeric(),end=numeric(),id=character(),variable=character())
  
  for (chosen_variable in soil_variables){
    soil_associations_combined %>% 
      filter(variable == chosen_variable) %>%
      filter(bayes_with_sv >= 20) -> bayes_20_hits
    
    for (n in 1:nrow(bayes_20_hits)){
      genes <- gff %>%
        filter(chr == bayes_20_hits$chr[n],
               (start - space_around_genes) <=bayes_20_hits$pos[n],
               (end + space_around_genes) >=bayes_20_hits$pos[n]) %>%
        mutate(variable = chosen_variable)
      bayes_20_genes <- rbind(bayes_20_genes, genes)
    }
  }
  for (chosen_variable in climate_variables){
    climate_associations_combined %>% 
      filter(variable == chosen_variable) %>%
      filter(bayes_with_sv >= 20) -> bayes_20_hits
    
    for (n in 1:nrow(bayes_20_hits)){
      genes <- gff %>%
        filter(chr == bayes_20_hits$chr[n],
               (start - space_around_genes) <=bayes_20_hits$pos[n],
               (end + space_around_genes) >=bayes_20_hits$pos[n]) %>%
        mutate(variable = chosen_variable)
      bayes_20_genes <- rbind(bayes_20_genes, genes)
    }
  }
  
  
  ###########
  #Save top regions for each variable
  ###########
  window_size <- 50
  climate_associations_windows <- climate_associations_combined %>%
    group_by(chr,variable) %>%
    mutate(counter = 1:n(),window=floor(counter/window_size)) %>%
    group_by(chr,variable,window) %>%
    
    summarize(start=min(pos),end=max(pos),median_bayes=median(bayes_with_sv)) %>%
    group_by(variable) %>%
    top_n(20, median_bayes)
  
  soil_associations_windows <- soil_associations_combined %>%
    group_by(chr,variable) %>%
    mutate(counter = 1:n(),window=floor(counter/window_size)) %>%
    group_by(chr,variable,window) %>%
    
    summarize(start=min(pos),end=max(pos),median_bayes=median(bayes_with_sv)) %>%
    group_by(variable) %>%
    top_n(20, median_bayes)
  
  write_tsv(climate_associations_windows, paste0("gea/Climate_GEA_",chosen_species,"_2019.tophits.txt"))
  write_tsv(soil_associations_windows, paste0("gea/Soil_GEA_",chosen_species,"_2019.tophits.txt"))
  
  ###########
  #Pick genes with bayes > 20 within windows
  ###########
  bayes_20_window_genes <- genes %>% head(0)
  for (chosen_variable in soil_variables){
    genes <- bayes_20_genes %>%
      filter(variable == chosen_variable)
    windows <- soil_associations_windows %>%
      filter(variable == chosen_variable)
    for (y in 1:nrow(windows)){
      chosen_chr <- windows$chr[y]
      chosen_start <- windows$start[y]
      chosen_end <- windows$end[y]
      left_overlap <- genes %>%
        filter(chr == chosen_chr) %>%
        filter(start <= chosen_start & end >= chosen_start)
      right_overlap <- genes %>%
        filter(chr == chosen_chr) %>%
        filter(start <= chosen_end & end >= chosen_end)
      middle_overlap <- genes %>%
        filter(chr == chosen_chr) %>%
        filter(start >= chosen_start & end <= chosen_end)
      bayes_20_window_genes <- rbind(bayes_20_window_genes, left_overlap,right_overlap, middle_overlap)
    }
  }
  for (chosen_variable in climate_variables){
    genes <- bayes_20_genes %>%
      filter(variable == chosen_variable)
    windows <- climate_associations_windows %>%
      filter(variable == chosen_variable)
    for (y in 1:nrow(windows)){
      chosen_chr <- windows$chr[y]
      chosen_start <- windows$start[y]
      chosen_end <- windows$end[y]
      left_overlap <- genes %>%
        filter(chr == chosen_chr) %>%
        filter(start <= chosen_start & end >= chosen_start)
      right_overlap <- genes %>%
        filter(chr == chosen_chr) %>%
        filter(start <= chosen_end & end >= chosen_end)
      middle_overlap <- genes %>%
        filter(chr == chosen_chr) %>%
        filter(start >= chosen_start & end <= chosen_end)
      bayes_20_window_genes <- rbind(bayes_20_window_genes, left_overlap,right_overlap, middle_overlap)
    }
  }
  write_tsv(unique(bayes_20_window_genes), paste0("gea/All_GEA_",chosen_species,"_2019.bayes20_top10windows.genes.txt"))
  next
  
  ###########
  #CLIMATE
  ###########
  
  #Save the Climate GEA plots without SVs
  pdf(paste0("gea/Climate_GEA_",chosen_species,"_2019.SVremoved.pdf"),height=6,width=12)
  for (i in 1:length(climate_variables)){
    chosen_variable = climate_variables[i]
    inversion_to_be_plotted <- climate_associations %>%
      filter(variable == chosen_variable) %>%
      inner_join(.,inv_locations)
    
    plot <- climate_associations_combined %>%
      filter(variable == chosen_variable) %>%
      filter(bayes_without_sv > 9) %>%
      inner_join(cum_location) %>%
      
      
      ggplot(.) +
      
      # Show all points
      geom_point( aes(x=poscum, y=bayes_without_sv,color=as.factor(chr)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
      geom_hline(yintercept = 10,linetype="dotted") +
      #geom_rect(data=ann_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0,ymax=5)) +
      
      # custom X axis:
      scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
      #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
      
      # Custom the theme:
      theme_bw() +
      theme( legend.position="none",
             panel.border = element_blank(),
             panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank()) +
      ggtitle(chosen_variable) +
      xlab("Chr")  +
      geom_segment(data=inversion_to_be_plotted %>%
                     filter(as.numeric(bayesfactor) > 10),
                   aes(x=startcum,xend=endcum,
                       y=as.numeric(bayesfactor),yend=as.numeric(bayesfactor)),
                   size=2,color="#9A2374",alpha=0.7) +
      ylab("dB")
    
    print(plot)
  }
  dev.off()
  
  #Save the Climate GEA plots with SVs
  pdf(paste0("gea/Climate_GEA_",chosen_species,"_2019.SVpresent.pdf"),height=6,width=12)
  for (i in 1:length(climate_variables)){
    chosen_variable = climate_variables[i]
    inversion_to_be_plotted <- climate_associations %>%
      filter(variable == chosen_variable) %>%
      inner_join(.,inv_locations)
    
    plot <- climate_associations_combined %>%
      filter(variable == chosen_variable) %>%
      filter(bayes_with_sv > 9) %>%
      inner_join(cum_location) %>%
      
      
      ggplot(.) +
      
      # Show all points
      geom_point( aes(x=poscum, y=bayes_with_sv,color=as.factor(chr)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
      geom_hline(yintercept = 10,linetype="dotted") +
      #geom_rect(data=ann_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0,ymax=5)) +
      
      # custom X axis:
      scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
      #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
      
      # Custom the theme:
      theme_bw() +
      theme( legend.position="none",
             panel.border = element_blank(),
             panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank()) +
      ggtitle(chosen_variable) +
      xlab("Chr")  +
      geom_segment(data=inversion_to_be_plotted %>%
                     filter(as.numeric(bayesfactor) > 10),
                   aes(x=startcum,xend=endcum,
                       y=as.numeric(bayesfactor),yend=as.numeric(bayesfactor)),
                   size=2,color="#9A2374",alpha=0.7) +
      ylab("dB")
    
    print(plot)
  }
  dev.off()
  
  #Climate compare with and without SVs
  pdf(paste0("gea/Climate_GEA_",chosen_species,"_2019.SVcompared.pdf"),height=6,width=12)
  for (i in 1:length(climate_variables)){
    chosen_variable = climate_variables[i]
    inversion_to_be_plotted <- climate_associations %>%
      filter(variable == chosen_variable) %>%
      inner_join(.,inv_locations)
    plot <- climate_associations_combined %>%
      filter(variable == chosen_variable) %>%
      mutate(bf_dif = bayes_without_sv - bayes_with_sv) %>%
      filter(abs(bf_dif) > 10) %>%
      inner_join(cum_location) %>%


      ggplot(.) +

      # Show all points
      geom_point( aes(x=poscum, y=bf_dif,color=as.factor(chr)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
      #geom_rect(data=ann_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0,ymax=5)) +

      # custom X axis:
      scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
      #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis

      # Custom the theme:
      theme_bw() +
      theme( legend.position="none",
             panel.border = element_blank(),
             panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank()) +
      ggtitle(chosen_variable) +
      xlab("Chr") +
      ylab("Difference in dB")

    print(plot)
   
  }
  dev.off()
  
  ###########
  #SOIL
  ###########
  #Save the Soil GEA plots without SVs
  pdf(paste0("gea/Soil_GEA_",chosen_species,"_2019.SVremoved.pdf"),height=6,width=12)
  for (i in 1:length(soil_variables)){
    chosen_variable = soil_variables[i]
    inversion_to_be_plotted <- soil_associations %>%
      filter(variable == chosen_variable) %>%
      inner_join(.,inv_locations)
    
    plot <- soil_associations_combined %>%
      filter(variable == chosen_variable) %>%
      filter(bayes_without_sv > 9) %>%
      inner_join(cum_location) %>%
      
      
      ggplot(.) +
      
      # Show all points
      geom_point( aes(x=poscum, y=bayes_without_sv,color=as.factor(chr)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
      geom_hline(yintercept = 10,linetype="dotted") +
      #geom_rect(data=ann_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0,ymax=5)) +
      
      # custom X axis:
      scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
      #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
      
      # Custom the theme:
      theme_bw() +
      theme( legend.position="none",
             panel.border = element_blank(),
             panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank()) +
      ggtitle(chosen_variable) +
      xlab("Chr")  +
      geom_segment(data=inversion_to_be_plotted %>%
                     filter(as.numeric(bayesfactor) > 10),
                   aes(x=startcum,xend=endcum,
                       y=as.numeric(bayesfactor),yend=as.numeric(bayesfactor)),
                   size=2,color="#9A2374",alpha=0.7) +
      ylab("dB")
    
    print(plot)
  }
  dev.off()
  
  #Save the soil GEA plots with SVs
  pdf(paste0("gea/Soil_GEA_",chosen_species,"_2019.SVpresent.pdf"),height=6,width=12)
  for (i in 1:length(soil_variables)){
    chosen_variable = soil_variables[i]
    inversion_to_be_plotted <- soil_associations %>%
      filter(variable == chosen_variable) %>%
      inner_join(.,inv_locations)
    
    plot <- soil_associations_combined %>%
      filter(variable == chosen_variable) %>%
      filter(bayes_with_sv > 9) %>%
      inner_join(cum_location) %>%
      
      
      ggplot(.) +
      
      # Show all points
      geom_point( aes(x=poscum, y=bayes_with_sv,color=as.factor(chr)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
      geom_hline(yintercept = 10,linetype="dotted") +
      #geom_rect(data=ann_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0,ymax=5)) +
      
      # custom X axis:
      scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
      #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
      
      # Custom the theme:
      theme_bw() +
      theme( legend.position="none",
             panel.border = element_blank(),
             panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank()) +
      ggtitle(chosen_variable) +
      xlab("Chr")  +
      geom_segment(data=inversion_to_be_plotted %>%
                     filter(as.numeric(bayesfactor) > 10),
                   aes(x=startcum,xend=endcum,
                       y=as.numeric(bayesfactor),yend=as.numeric(bayesfactor)),
                   size=2,color="#9A2374",alpha=0.7) +
      ylab("dB")
    
    print(plot)
  }
  dev.off()
  
  #Climate compare with and without SVs
  pdf(paste0("gea/Soil_GEA_",chosen_species,"_2019.SVcompared.pdf"),height=6,width=12)
  for (i in 1:length(soil_variables)){
    chosen_variable = soil_variables[i]
    inversion_to_be_plotted <- soil_associations %>%
      filter(variable == chosen_variable) %>%
      inner_join(.,inv_locations)
    plot <- soil_associations_combined %>%
      filter(variable == chosen_variable) %>%
      mutate(bf_dif = bayes_without_sv - bayes_with_sv) %>%
      filter(abs(bf_dif) > 10) %>%
      inner_join(cum_location) %>%
      
      
      ggplot(.) +
      
      # Show all points
      geom_point( aes(x=poscum, y=bf_dif,color=as.factor(chr)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
      #geom_rect(data=ann_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0,ymax=5)) +
      
      # custom X axis:
      scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
      #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
      
      # Custom the theme:
      theme_bw() +
      theme( legend.position="none",
             panel.border = element_blank(),
             panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank()) +
      ggtitle(chosen_variable) +
      xlab("Chr") +
      ylab("Difference in dB")
    
    print(plot)
    
  }
  dev.off()
  
}
##########
###TRYING OUT WAYS TO PICK OUTLIERS
##########
climate_associations_with_sv %>%
  gather(.,variable,bayes_with_sv,latitude_e:RH_e)  %>%
  filter(variable != "PAS_e") -> climate_associations_combined

window_size <- 50 #Number of SNPs in window
climate_associations_windows <- climate_associations_combined %>%
  group_by(chr,variable) %>%
  mutate(counter = 1:n(),window=floor(counter/window_size)) %>%
  group_by(chr,variable,window) %>%
  
  summarize(start=min(pos),end=max(pos),median_bayes=median(bayes_with_sv)) %>%
  group_by(variable) %>%
  top_n(10, median_bayes)

climate_associations_windows %>%
  group_by(variable) %>%
  top_n(10, median_bayes  )
climate_associations_windows_cum <- chr_lengths %>% 
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(climate_associations_windows, by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate( startcum=start+tot,endcum=end+tot)
axisdf = climate_associations_windows_cum %>% group_by(chr) %>% summarize(center=( max(startcum) + min(startcum) ) / 2 )


climate_associations_windows_cum %>%
  filter(variable == "elevation_e") %>%
  ggplot(.) +
  
  # Show all points
  geom_point( aes(x=startcum, y=median_bayes,color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  #geom_rect(data=ann_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0,ymax=5)) +
  
  # custom X axis:
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( legend.position="none",
         panel.border = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank()) +
  xlab("Chr") +
  ylab("Difference in dB")

climate_associations_combined %>%
  group_by(variable) %>%
  top_n(100, bayes_with_sv  ) -> climate_associations_combined_top100

climate_associations_combined_top100 %>%
  group_by(variable,chr) %>%
  mutate(pos_dif_rev = pos - lag(pos),
         pos_dif_for = lead(pos) - pos) %>%
  mutate(pos_dif = case_when(is.na(pos_dif_rev) ~ pos_dif_for,
                             is.na(pos_dif_for) ~ pos_dif_rev,
                             pos_dif_rev < pos_dif_for ~ pos_dif_rev,
                             pos_dif_rev >= pos_dif_for ~ pos_dif_for)) %>%
  filter(pos_dif < 100000) %>% View()
  


