#Printing pleotropy of GWAS and GEA
library(tidyverse)
library(ggrastr)

pvalue_cutoff <- 0.001
min_freq <- 0.03
directory <- "/media/owens/Copper/wild_gwas/gwas"
species_abbreviation_list <- c("annuus","argophyllus","petpet","petfal")
species_capital_list <-c("Annuus", "Argophyllus","Petiolaris","Petiolaris")
tag_list <- c("gwas", "gwas","petpet","petfal")
trait_types <- read_tsv("gwas/all_traits_category.txt",col_names = c("trait","trait_type"))
######Compress the bigger table down to a smaller version for the main text

all_gwas_inv <- read_tsv(paste0(directory,"/","annuus","/results/",
                                "Ha412HO_inv.v3.pcasites.","annuus",
                                ".shuffled.",
                                "DTF",".ps.gz"),
                         col_names=c("id","beta","sd","pvalue")) %>% mutate(trait = "DTF",
                                                                            species = "annuus") %>%
  head(0)
all_gea_inv <- tibble(type=character(),sig_count=numeric(),id=character(),species=character())

for (n in 1:4){
  
  
  species_abbreviation <- species_abbreviation_list[n]
  species_capital <- species_capital_list[n]
  tag <- tag_list[n]
  
  ##Get list of all traits
  traits <- tibble(trait=list.files(paste0("gwas/",species_abbreviation),include.dirs = F) ) %>%
    filter(!grepl("znorm",trait)) %>%
    separate(trait,c("trait","spacer"),"\\.") %>%
    select(-spacer)
  
  all_species_inv <- read_tsv(paste0(directory,"/","annuus","/results/",
                                     "Ha412HO_inv.v3.pcasites.","annuus",
                                     ".shuffled.",
                                     "DTF",".ps.gz"),
                              col_names=c("id","beta","sd","pvalue")) %>% mutate(trait = traits$trait[1],
                                                                                 species = species_abbreviation) %>%
    head(0)
  
  #Pull out associations for inversions
  
  all_invs <- read_tsv(paste0(directory,"/",species_abbreviation,"/results/",
                              "Ha412HO_inv.v3.pcasites.",species_abbreviation,
                              ".shuffled.",
                              traits$trait[1],".ps.gz"),
                       col_names=c("id","beta","sd","pvalue")) %>%
    select(id) %>% separate(id,c("chr","pos"),"_") %>%
    mutate(count=0)
  
  for (x in 1:nrow(traits)){
    gwas_inv <- read_tsv(paste0(directory,"/",species_abbreviation,"/results/",
                                "Ha412HO_inv.v3.pcasites.",species_abbreviation,
                                ".shuffled.",
                                traits$trait[x],".ps.gz"),
                         col_names=c("id","beta","sd","pvalue")) %>%
      mutate(trait = traits$trait[x],
             species = species_abbreviation)
    all_species_inv <- rbind(all_species_inv,gwas_inv)
    
    
    
  }
  all_species_inv %>%  mutate(id = gsub("_", ".0",id)) -> all_species_inv
  
  ##Remove inversions below a minimum allele frequency
  if (species_capital == "Petiolaris"){
    inv_frequencies <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.petinvfreq.txt") %>%
      filter(tolower(species) == species_abbreviation) %>%
      mutate(chr = gsub("Ha412HOChr","",chr),
             pos=gsub("syn","",mds)) %>%
      select(-mds,-species) %>%
      mutate(id = paste0(chr,"_",pos)) %>%
      mutate(id = gsub("_",".0",id))
    
    inner_join(all_species_inv, inv_frequencies) %>%
      filter(freq > min_freq,freq < (1-min_freq)) %>%
      select(-freq,-chr,-pos) -> all_species_inv
    
  }
  all_gwas_inv <- rbind(all_gwas_inv,all_species_inv)

    
  
  #Now load in GEA
  gea_traits <- read_tsv("gea/environment_variables.txt") 
  gea <- read_tsv(paste0("/media/owens/Copper/wild_gwas/env_associations/climate/sv_HA412/",species_abbreviation,
                           "/varout_",species_abbreviation,"_HA412_remapped_Baypass_noinv_matrix_25_vars.tab")) %>% 
    separate(inversion_ID, c("chr","mds")) %>%
    mutate(chr_n = gsub("Ha412HOChr","",chr)) %>% 
    select(-chr) %>%
    gather(variable, bayesfactor, MAT_e:RH_e) %>%
    inner_join(gea_traits) %>%
    group_by(chr_n,mds,type) %>%
    mutate(bayesfactor = as.numeric(bayesfactor)) %>%
    summarize(sig_count = sum(bayesfactor > 10),
              max_sig = max(bayesfactor)) %>%
    mutate(any_sig = case_when(sig_count > 0 ~ "significant",
                               TRUE ~ "nonsignificant")) %>%
    ungroup() %>%
    mutate(id = paste0(chr_n,"_",mds)) %>%
    mutate(id = gsub("syn","",id)) %>%
    mutate(id = gsub("pos","",id)) %>%
    mutate(id = gsub("neg","",id)) %>%
    mutate(id = gsub("_", ".0",id)) %>%
    mutate(species = species_abbreviation) %>%
    select(type,sig_count,id,species)
  
  soil_associations <-read_tsv(paste0("/media/owens/Copper/wild_gwas/env_associations/soil/sv_HA412/baypass_inversions_soil_data_bf.txt")) %>%
    filter(species == species_abbreviation) %>%
    mutate(full_id = paste0(chr,"_",mds)) %>%
    mutate(inversion_ID = full_id) %>%
    mutate(inversion_ID = gsub("Ha412HOChr","",inversion_ID)) %>%
    mutate(inversion_ID = gsub("pos","",inversion_ID)) %>%
    mutate(inversion_ID = gsub("neg","",inversion_ID)) %>%
    mutate(inversion_ID = gsub("syn","",inversion_ID)) %>%
    mutate(inversion_ID = gsub("_",".0",inversion_ID)) %>%
    select(-ID,-species,-mds) %>%
    gather(variable,bayesfactor,colnames(.)[2]:colnames(.)[ncol(.)-2]) %>%
    select(-chr,-full_id) %>%
    inner_join(gea_traits) %>%
    group_by(inversion_ID) %>%
    mutate(bayesfactor = as.numeric(bayesfactor)) %>%
    summarize(sig_count = sum(bayesfactor > 10),
              max_sig = max(bayesfactor)) %>%
    mutate(any_sig = case_when(sig_count > 0 ~ "significant",
                               TRUE ~ "nonsignificant")) %>%
    ungroup() %>%
    rename(id=inversion_ID) %>%
    mutate(type = "soil") %>%
    mutate(species = species_abbreviation) %>%
    select(type,sig_count,id,species)
  
  gea <- rbind(gea, soil_associations)
  
  if (species_capital == "Petiolaris"){
    inv_frequencies <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.petinvfreq.txt") %>%
      filter(tolower(species) == species_abbreviation) %>%
      mutate(chr = gsub("Ha412HOChr","",chr),
             pos=gsub("syn","",mds)) %>%
      select(-mds,-species) %>%
      mutate(id = paste0(chr,"_",pos)) %>%
      mutate(id = gsub("_",".0",id))
    inner_join(gea, inv_frequencies) %>%
      filter(freq > min_freq,freq < (1-min_freq)) %>%
      select(-freq,-chr,-pos) -> gea
    
  }
  all_gea_inv <- rbind(all_gea_inv, gea)
}
all_gwas_inv %>%
  mutate(new_id = paste0(species,"_",id)) %>%
  select(species,new_id) %>%
  unique() %>%
  arrange(species) %>%
  mutate(count=1:nrow(.)) %>%
  group_by(species) %>%
  summarize(start=min(count),end=max(count)) -> background


pdf("Ha412HO_inv.v3.pcasites.gwas.gea.typesummary.pdf",height=3,width=8)
all_gwas_inv %>%
  inner_join(trait_types) %>%
  mutate(significant = case_when(pvalue < pvalue_cutoff ~ "1",
                                 TRUE ~ "0")) %>%
  group_by(id,trait_type,species) %>%
  summarise(any_sig = sum(as.numeric(significant))) %>%
  mutate(any_sig = case_when(any_sig >= 1 ~ 1,
                             TRUE ~ 0)) %>%
  filter(trait_type != "Pigmentation") %>%
  ggplot(.,) +
  geom_tile(aes(x=paste0(species,id),y=fct_rev(trait_type),alpha=any_sig)) +
  #facet_wrap(~species,scales="free_x") +
  theme_minimal() +
  ylab("Associated trait categories") + xlab("SV") +
  annotate("rect",xmin=background %>% filter(species == "annuus") %>% pull(start)-0.5,
           xmax=background %>% filter(species == "annuus") %>% pull(end)+0.5,
           ymin=0.5,ymax=6.5,fill="#FFB31C",alpha=0.6) +
  annotate("rect",xmin=background %>% filter(species == "argophyllus") %>% pull(start)-0.5,
           xmax=background %>% filter(species == "argophyllus") %>% pull(end)+0.5,
           ymin=0.5,ymax=6.5,fill="#3DA047",alpha=0.6) +
  annotate("rect",xmin=background %>% filter(species == "petfal") %>% pull(start)-0.5,
           xmax=background %>% filter(species == "petfal") %>% pull(end)+0.5,
           ymin=0.5,ymax=6.5,fill="#447499",alpha=0.6) +
  annotate("rect",xmin=background %>% filter(species == "petpet") %>% pull(start)-0.5,
           xmax=background %>% filter(species == "petpet") %>% pull(end)+0.5,
           ymin=0.5,ymax=6.5,fill="#4BB7B7",alpha=0.6) +
  geom_tile(data=all_gwas_inv %>%
              inner_join(trait_types) %>%
              mutate(significant = case_when(pvalue < pvalue_cutoff ~ "1",
                                             TRUE ~ "0")) %>%
              group_by(id,trait_type,species) %>%
              summarise(any_sig = sum(as.numeric(significant))) %>%
              mutate(any_sig = case_when(any_sig >= 1 ~ 1,
                                         TRUE ~ 0)) %>%
              filter(trait_type != "Pigmentation") %>%
              filter(any_sig > 0),
            aes(x=paste0(species,id),y=fct_rev(trait_type)),fill="black") +
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(
        labels=all_gwas_inv %>% group_by(species, id) %>% summarize(count=n()) %>% mutate(label = paste0(substr(species,0,3),id)) %>% pull(label)) +
  annotate("rect",xmin=background %>% filter(species == "argophyllus") %>% pull(start)-0.5,
           xmax=background %>% filter(species == "argophyllus") %>% pull(end)+0.5,
           ymin=0.5,ymax=1.5,fill="white",alpha=0.6) +
  annotate("rect",xmin=background %>% filter(species == "argophyllus") %>% pull(start)-0.5,
           xmax=background %>% filter(species == "argophyllus") %>% pull(end)+0.5,
           ymin=2.5,ymax=3.5,fill="white",alpha=0.6) +
  annotate("rect",xmin=background %>% filter(species == "argophyllus") %>% pull(start)-0.5,
           xmax=background %>% filter(species == "argophyllus") %>% pull(end)+0.5,
           ymin=5.5,ymax=6.5,fill="white",alpha=0.6) +
  geom_vline(xintercept=background$end+0.5,color="white")


gea_plotting_data <- all_gea_inv %>%
  mutate(any_sig = case_when(sig_count >= 1 ~ 1,
                             TRUE ~ 0)) %>%
  rename(trait_type = type)   %>% 
  mutate(trait_type = case_when(trait_type == "soil" ~ "Soil",
          trait_type == "Combined" ~ "T+P",
          TRUE ~ trait_type))
gea_plotting_data$trait_type <- gea_plotting_data$trait_type %>%
  fct_relevel(., c("Soil", "Precipitation","Temperature","P+T"))
  
gea_plotting_data %>%

  ggplot(.,) +
  geom_tile(aes(x=paste0(species,id),y=fct_rev(trait_type),alpha=any_sig)) +
  #facet_wrap(~species,scales="free_x") +
  theme_minimal() +
  ylab("Associated climate categories") + xlab("SV") +
  annotate("rect",xmin=background %>% filter(species == "annuus") %>% pull(start)-0.5,
           xmax=background %>% filter(species == "annuus") %>% pull(end)+0.5,
           ymin=0.5,ymax=4.5,fill="#FFB31C",alpha=0.6) +
  annotate("rect",xmin=background %>% filter(species == "argophyllus") %>% pull(start)-0.5,
           xmax=background %>% filter(species == "argophyllus") %>% pull(end)+0.5,
           ymin=0.5,ymax=4.5,fill="#3DA047",alpha=0.6) +
  annotate("rect",xmin=background %>% filter(species == "petfal") %>% pull(start)-0.5,
           xmax=background %>% filter(species == "petfal") %>% pull(end)+0.5,
           ymin=0.5,ymax=4.5,fill="#447499",alpha=0.6) +
  annotate("rect",xmin=background %>% filter(species == "petpet") %>% pull(start)-0.5,
           xmax=background %>% filter(species == "petpet") %>% pull(end)+0.5,
           ymin=0.5,ymax=4.5,fill="#4BB7B7",alpha=0.6) +
  geom_tile(data=gea_plotting_data %>%
              mutate(any_sig = case_when(sig_count >= 1 ~ 1,
                                         TRUE ~ 0)) %>%
              filter(any_sig > 0),
            aes(x=paste0(species,id),y=fct_rev(trait_type)),fill="black") +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(labels=all_gwas_inv %>% group_by(species, id) %>% summarize(count=n()) %>% mutate(label = paste0(substr(species,0,3),id)) %>% pull(label)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(xintercept=background$end+0.5,color="white")
dev.off()
