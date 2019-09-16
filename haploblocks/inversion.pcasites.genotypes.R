#Plotting fst outliers, but actual calls.

library(tidyverse)
library(ggthemes)
library(forcats)
library(ggrastr)
library(gridExtra)
library(scatterpie)

#Load in all inversion genotypes
pca_genotypes <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.genotypes.v1.txt")
folder <- "MDS_outliers"
chosen_species <- "petiolaris"
chosen_species_file <- "petiolaris_syn"
chosen_species_abbreviation <- c("PetPet","PetFal")
filtering <- "pcasites"
prefix <- "Ha412HO_inv.v3"
base_directory <- "/media/owens/Copper/wild_gwas/"
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T)
pop_loc <- read_tsv("pop_loc_allnum.txt") %>% rename(population = pop)
pca_genotypes <- pca_genotypes %>% filter(species == chosen_species)
inversion_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt") %>%
  rename(chr_n = chr) %>%
  mutate(chr = paste0("Ha412HOChr",sprintf("%02d", chr_n)),
         mds = paste0(direction,mds))
inversions <- inversion_list %>% filter(spe == chosen_species_file) %>%
  select(chr,mds) %>% unique()


for (n in (1:nrow(inversions))){
  chosen_chr <- inversions[n,1] %>% pull()
  chosen_mds <- inversions[n,2] %>% pull()
  
  if (chosen_species == "annuus"){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(25, 50)
    long_range <- c(-125,-93)
    pie_size <- 0.4
  }else if(chosen_species == "argophyllus" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- subset(states, region %in% c("texas"))
    lat_range <- c(25, 30)
    long_range <- c(-99,-96)
    pie_size = 0.1
  }else if(chosen_species == "petiolaris" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(30, 43)
    long_range <- c(-115,-95)
    pie_size <- 0.4
  }
  
  
  
  fst_calls <- read_tsv(paste(base_directory,chosen_species,"/",prefix,".",chosen_species_file,".", chosen_chr,".",chosen_mds,".",
                filtering,".genotypes.txt.gz",sep=""),guess_max= 100000) %>%
    inner_join(., labels %>% dplyr::rename(sample=name)) %>%
    filter(genotype == "00" | genotype == "01" | genotype == "11")
  
  min_depth_genotyping <- 5
  min_depth_visualizing <- 2
  fst_calls %>%
    filter(depth >= min_depth_genotyping) %>%
    group_by(sample, species,population, is_wild, genotype) %>%
    count() %>%
    filter(!is.na(genotype)) %>%
    spread(genotype, n,fill=0) %>%
    mutate(total =  (`00` + `01` + `11`)*2,
           perc_1 = ((`11`*2) + `01`)/total,
           het = (`01`*2)/total,
           left_formula = ((-(2/3)*perc_1) + (2/3)),
           right_formula = (((2/3)*perc_1))) %>%
    mutate(fst_call = case_when(perc_1 < 0.5  & left_formula >= het ~ 0,
                                perc_1 < 0.5  & left_formula < het ~ 1,
                                perc_1 >= 0.5  & right_formula >= het ~ 2,
                                perc_1 >= 0.5  & right_formula < het ~ 1)) -> fst_samples
  

  if (exists("pca_genotypes")){
  pca_genotypes %>% 
    filter(chr == chosen_chr,mds== chosen_mds) %>%
    select(sample,cluster) %>%
    full_join(fst_samples) %>% 
    select(sample, cluster, total, perc_1, het, fst_call) %>%
    rename(cluster_genotype = cluster, triangle_genotype = fst_call, pca_sites = total) -> output_genotypes
  }else{
    fst_samples %>%
      mutate(cluster = "NA") %>%
      select(sample, cluster, total, perc_1, het, fst_call) %>% 
      rename(cluster_genotype = cluster, triangle_genotype = fst_call, pca_sites = total) %>%
      ungroup() %>%
      select(sample, cluster_genotype, pca_sites, perc_1, het, triangle_genotype) -> output_genotypes
  }
  write_tsv(output_genotypes, paste(folder,"/Ha412HO/",chosen_species_file,"/",
                                    prefix,".",filtering,".",chosen_chr,".",chosen_mds,".genotypes.txt",sep=""))


  
  pdf(paste(folder,"/Ha412HO/",chosen_species_file,"/",prefix,".",filtering,".",chosen_chr,".",chosen_mds,".genotypes.pdf",sep=""),
      width=12,height=12)
  

    if (exists("pca_genotypes")){
      print(
      pca_genotypes %>% filter(chr == chosen_chr,mds== chosen_mds) %>%
        select(sample,cluster) %>%
        inner_join(.,fst_samples) %>%
        mutate(agreement = case_when(fst_call == cluster ~ "Concordant",
                                     TRUE ~ "Discordant")) %>%
        ggplot(.,aes(species, fill=agreement)) + geom_bar(position="stack") +
        theme_bw() +
        scale_fill_brewer(palette = "Set1")
      )
    }

  print(
  fst_samples %>%
    filter(total > 10) %>%
    ggplot(.,aes(x=perc_1,y=het)) +
    geom_segment(aes(x=0,xend=0.5,y=0,yend=1),color="grey") +
    geom_segment(aes(x=0.5,xend=1,y=1,yend=0),color="grey") +
    geom_segment(aes(x=0,xend=1,y=0,yend=0),color="grey") +
    geom_point(alpha=0.4) +
    theme_few() +
    ylab("Heterozygosity") +
    xlab("Proportion `1` allele") +
    facet_wrap(~species) +
    ggtitle(paste(chosen_species,chosen_chr,chosen_mds))
  )
  
  for (x in 1:length(chosen_species_abbreviation)){
    print(
      fst_samples %>%
        filter(total > 10) %>%
        filter(species == chosen_species_abbreviation[x]) %>%
        ggplot(.,aes(x=perc_1,y=het)) +
        geom_segment(aes(x=0,xend=0.5,y=0,yend=1),color="grey") +
        geom_segment(aes(x=0.5,xend=1,y=1,yend=0),color="grey") +
        geom_segment(aes(x=0,xend=1,y=0,yend=0),color="grey") +
        geom_point(alpha=0.4,aes(color=as.factor(fst_call))) +
        geom_segment(aes(x=0.25,xend=0.5,y=0.5,yend=(1/3)),linetype="dotted") +
        geom_segment(aes(x=0.5,xend=0.75,y=(1/3),yend=0.5),linetype="dotted") +
        geom_segment(aes(x=0.5,xend=0.5,y=0,yend=(1/3)),linetype="dotted") +
        scale_color_brewer(palette = "Set1",name="Genotype") +
        theme_few() +
        ylab("Heterozygosity") +
        xlab("Proportion `1` allele") +
        ggtitle(paste(chosen_species,chosen_chr,chosen_mds)) +
        facet_wrap(~species)
    )
    if (chosen_species == "annuus"){
      x <- 1
      print(
        fst_samples %>%
          filter(total > 10) %>%
          filter(species == chosen_species_abbreviation[x]) %>%
          filter(is_wild != "wild") %>%
          ggplot(.,aes(x=perc_1,y=het)) +
          geom_segment(aes(x=0,xend=0.5,y=0,yend=1),color="grey") +
          geom_segment(aes(x=0.5,xend=1,y=1,yend=0),color="grey") +
          geom_segment(aes(x=0,xend=1,y=0,yend=0),color="grey") +
          geom_point(alpha=0.4,aes(color=as.factor(fst_call))) +
          geom_segment(aes(x=0.25,xend=0.5,y=0.5,yend=(1/3)),linetype="dotted") +
          geom_segment(aes(x=0.5,xend=0.75,y=(1/3),yend=0.5),linetype="dotted") +
          geom_segment(aes(x=0.5,xend=0.5,y=0,yend=(1/3)),linetype="dotted") +
          scale_color_brewer(palette = "Set1",name="Genotype") +
          theme_few() +
          ylab("Heterozygosity") +
          xlab("Proportion `1` allele") +
          ggtitle(paste(chosen_species,chosen_chr,chosen_mds)) +
          facet_wrap(~is_wild)
      )
      
    }
    fst_calls %>%
      filter(depth >= min_depth_visualizing) %>%
      group_by(sample) %>%
      filter(species == chosen_species_abbreviation[x],is_wild == "wild") %>%
      mutate(count = n()) %>%
      ungroup() %>%
      mutate(max_count = max(count)) %>% 
      filter(count > max_count/10) %>%
      inner_join(., fst_samples) %>% select(pos) %>% unique() %>% pull() -> positions
    n_labels <- ceiling(length(positions)/40)
    x_labels <- sort(positions)[seq(0, length(positions), by= n_labels)]
    
    print(
      fst_calls %>%
        filter(depth >= min_depth_visualizing) %>%
        group_by(sample) %>%
        filter(species == chosen_species_abbreviation[x],is_wild == "wild") %>%
        mutate(count = n()) %>%
        ungroup() %>%
        mutate(max_count = max(count)) %>% 
        filter(count > max_count/10) %>%
        inner_join(., fst_samples) %>% 
        ggplot(.,aes(y=fct_reorder(sample, perc_1),x=as.factor(pos),fill=as.factor(genotype))) + 
        geom_tile() + 
        theme_few() + 
        scale_fill_brewer(palette = "Set1",name="Site Genotype") +
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        ylab("Sample") + xlab("Position") +
        ggtitle(paste("min depth =", min_depth_visualizing)) +
        scale_x_discrete(breaks=x_labels) +
        theme(axis.text.x=element_text(angle=60, hjust=1)) + 
        facet_wrap(~species)
    )
    if (chosen_species == "annuus"){
      fst_calls %>%
        filter(depth >= min_depth_visualizing) %>%
        group_by(sample) %>%
        filter(species == chosen_species_abbreviation[x],is_wild != "wild") %>%
        mutate(count = n()) %>%
        ungroup() %>%
        mutate(max_count = max(count)) %>% 
        filter(count > max_count/10) %>%
        inner_join(., fst_samples) %>% select(pos) %>% unique() %>% pull() -> positions
      n_labels <- ceiling(length(positions)/40)
      x_labels <- sort(positions)[seq(0, length(positions), by= n_labels)]
      
      print(
        fst_calls %>%
          filter(depth >= min_depth_visualizing) %>%
          group_by(sample) %>%
          filter(species == chosen_species_abbreviation[x],is_wild != "wild") %>%
          mutate(count = n()) %>%
          ungroup() %>%
          mutate(max_count = max(count)) %>% 
          filter(count > max_count/10) %>%
          inner_join(., fst_samples) %>% 
          ggplot(.,aes(y=fct_reorder(sample, perc_1),x=as.factor(pos),fill=as.factor(genotype))) + 
          geom_tile() + 
          theme_few() + 
          scale_fill_brewer(palette = "Set1",name="Site Genotype") +
          theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank()) +
          ylab("Sample") + xlab("Position") +
          ggtitle(paste("min depth =", min_depth_visualizing)) +
          scale_x_discrete(breaks=x_labels) +
          theme(axis.text.x=element_text(angle=60, hjust=1)) + 
          facet_wrap(~is_wild,nrow=2,scales = "free_y")
      )
    }
  }
  print(
  ggplot(target_state, aes(long, lat)) +
    geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
    coord_quickmap() +
    geom_scatterpie(data=inner_join(output_genotypes, info) %>% 
                      inner_join(., pop_loc) %>%
                      filter(species %in% chosen_species_abbreviation) %>%
                      group_by(population, lat, long, triangle_genotype) %>% 
                      tally() %>%
                      spread(., triangle_genotype, n,fill=0),
                    aes(x=long, y=lat, r=pie_size), 
                    cols=c("0","1","2"), color=NA, alpha=.8) +
    scale_fill_brewer(name="Cluster",palette = "Set1") +theme_bw() +
    xlab("Longitude") + ylab("Latitude") +
    scale_x_continuous(limits = long_range, expand = c(0, 0)) +
    scale_y_continuous(limits = lat_range, expand = c(0, 0))
  )
  
  dev.off()
}




###Niveus chr5 weirdness
pdf(paste(folder,"/Ha412HO/","niveus","/",prefix,".",filtering,".",chosen_chr,".",chosen_mds,".allspecies.genotypes.pdf",sep=""),
    width=12,height=20)
fst_calls %>%
  filter(depth >= min_depth_visualizing) %>%
  group_by(sample) %>%
  #filter(species == chosen_species_abbreviation[x],is_wild != "wild") %>%
  mutate(count = n()) %>%
  ungroup() %>%
  mutate(max_count = max(count)) %>% 
  filter(count > max_count/10) %>%
  inner_join(., fst_samples) %>% 
  ggplot(.,aes(y=fct_reorder(sample, perc_1),x=as.factor(pos),fill=as.factor(genotype))) + 
  geom_tile() + 
  theme_few() + 
  scale_fill_brewer(palette = "Set1",name="Site Genotype") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Sample") + xlab("Position") +
  ggtitle(paste("min depth =", min_depth_visualizing)) +
  scale_x_discrete(breaks=x_labels) +
  theme(axis.text.x=element_text(angle=60, hjust=1)) + 
  facet_wrap(~species,ncol=2,scales = "free_y")
dev.off()


####

fst_calls %>%
  filter(depth >= min_depth) %>%
  group_by(sample, species,population, is_wild, genotype) %>%
  count() %>%
  filter(!is.na(genotype)) %>%
  spread(genotype, n,fill=0) %>%
  mutate(total =  (`00` + `01` + `11`)*2,
         perc_1 = ((`11`*2) + `01`)/total,
         het = (`01`*2)/total) %>%
  mutate(fst_call = case_when(perc_1 < 0.1 ~ 0,
                              perc_1 > 0.9 ~ 2,
                              TRUE ~ 1)) -> fst_samples




pdf(paste(prefix, ".mdsoutlierfst.calls.pdf",sep=""),width=8,height=6)
for (mds_chosen in sort(unique(fst_samples$mds_coord))){
  if(fst_samples %>%filter(total > 40, mds_coord == mds_chosen) %>% nrow() == 0){next}
  fst_samples %>%
    filter(total > 40) %>%
    ggplot(.,aes(x=perc_1,y=het)) +
    geom_segment(aes(x=0,xend=0.5,y=0,yend=1),color="grey") +
    geom_segment(aes(x=0.5,xend=1,y=1,yend=0),color="grey") +
    geom_segment(aes(x=0,xend=1,y=0,yend=0),color="grey") +
    geom_point() +
    theme_few() +
    ylab("Heterozygosity") +
    xlab("Proportion `1` allele") +
    facet_wrap(~species)-> p1
  
  fst_samples %>%
    #filter(total > 40, mds_coord == mds_chosen) %>%
    ggplot(.,aes(perc_1)) +
    geom_histogram() +
    theme_few() +
    ylab("Count") +
    xlab("Proportion `1` allele") +
    facet_wrap(~species,scales="free_y")  ->p2
  print(  
    grid.arrange(p1, p2, nrow = 2,top=mds_chosen)
  )
}
dev.off()


pdf(paste(prefix, ".mdsoutlierfst.calls.justann.pdf",sep=""),width=8,height=6)
for (mds_chosen in sort(unique(fst_samples$mds_coord))){
  if(fst_samples %>%filter(total > 20, mds_coord == mds_chosen) %>% nrow() == 0){next}
  fst_samples %>%
    filter(species == "Ann") %>%
    filter(total > 20) %>%
    ggplot(.,aes(x=perc_1,y=het)) +
    geom_segment(aes(x=0,xend=0.5,y=0,yend=1),color="grey") +
    geom_segment(aes(x=0.5,xend=1,y=1,yend=0),color="grey") +
    geom_segment(aes(x=0,xend=1,y=0,yend=0),color="grey") +
    geom_point() +
    theme_few() +
    ylab("Heterozygosity") +
    xlab("Proportion `1` allele") +
    facet_wrap(~is_wild)-> p1
  
  fst_samples %>%
    filter(species == "Ann") %>%
    filter(total > 20, mds_coord == mds_chosen) %>%
    ggplot(.,aes(perc_1)) +
    geom_histogram() +
    theme_few() +
    ylab("Count") +
    xlab("Proportion `1` allele") +
    geom_vline(data=fst_samples %>%
                 filter(name == "SAM053", mds_coord == mds_chosen),
               aes(xintercept=perc_1),color="blue",alpha=0.8) +
    geom_vline(data=fst_samples %>%
                 filter(name == "SAM261", mds_coord == mds_chosen),
               aes(xintercept=perc_1),color="red",alpha=0.8) +
    facet_wrap(~is_wild,scales="free_y") ->p2
  print(  
    grid.arrange(p1, p2, nrow = 2,top=mds_chosen)
  )
}
dev.off()


pdf(paste(prefix, ".mdsoutlierfst.sitecalls.pdf",sep=""),width=5,height=10)
for (mds_chosen in sort(unique(fst_calls$mds_coord))){
  if (fst_calls %>%
      filter(depth >= min_depth) %>%
      filter(mds_coord == mds_chosen,!is.na(chr)) %>% nrow() == 0){next()}
  print(
    fst_calls %>%
      filter(depth >= min_depth) %>%
      group_by(sample) %>%
      filter(species == "Ann",is_wild == "wild") %>%
      mutate(count = n()) %>%
      ungroup() %>%
      mutate(max_count = max(count)) %>% 
      filter(count > max_count/10) %>%
      inner_join(., fst_samples) %>% 
      ggplot(.,aes(y=fct_reorder(sample, perc_1),x=as.factor(pos),fill=as.factor(genotype))) + 
      geom_tile() + 
      theme_few() + 
      scale_fill_brewer(palette = "Set1",name="Site Genotype") +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      ylab("Sample") + xlab("Position") +
      ggtitle(paste("min depth =", min_depth)) 
      #facet_grid(population~.,scales="free_y")
  )
}
dev.off()

pdf(paste(prefix, ".mdsoutlierfst.sitecalls.sam.pdf",sep=""),width=5,height=30)
fst_samples %>%
  filter(is_wild == "cultivar") %>%
  filter(total > 20) %>% 
  mutate(call = case_when(perc_1 < 0.15 ~ "0",
                          perc_1 > 0.85 ~ "2",
                          TRUE ~ "1"))%>%
  ggplot(.) +
  geom_tile(aes(y=name,x=mds_coord,fill=call)) +
  scale_fill_brewer(palette = "Set1",name="Inversion call") +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=5))
dev.off()
