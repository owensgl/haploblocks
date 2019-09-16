library(tidyverse)
library(ggrepel)

species_list <- c("annuus","argophyllus","petfal","petpet")
capital_species_list <-  c("Annuus","Argophyllus","Petiolaris","Petiolaris")
chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")
min_freq <- 0.03

inversion_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt") %>%
  mutate(inversion_ID = paste0("Ha412HOChr",sprintf("%02d",chr),"_",direction,mds)) %>%
  select(inversion_ID, sv_name,species)

all_associations <- tibble(sv_name=character(),variable=character(),bayesfactor=numeric(),sv_species_name=character())
#pdf("gea/Ha412HO.climategea.sv.baypass.pdf",height=15,width=12)
for (n in 1:4){
  chosen_species <- species_list[n]
  chosen_capital_species <- capital_species_list[n]
  
  env_associations <-  read_tsv(paste0("/media/owens/Copper/wild_gwas/env_associations/climate/sv_HA412/",chosen_species,
                                       "/varout_",chosen_species,"_HA412_remapped_Baypass_noinv_matrix_25_vars.tab")) %>% 
    inner_join(inversion_list %>% filter(species == chosen_capital_species)) %>% 
    select(-species,-inversion_ID) %>%
    gather(variable,bayesfactor,-sv_name) %>%
    mutate(bayesfactor = as.numeric(bayesfactor))
  
  soil_associations <-read_tsv(paste0("/media/owens/Copper/wild_gwas/env_associations/soil/sv_HA412/baypass_inversions_soil_data_bf.txt")) %>%
    filter(species == chosen_species) %>% 
    mutate(inversion_ID = paste0(chr,"_",mds)) %>% 
    select(-chr,-ID,-species,-mds) %>%
    
    inner_join(inversion_list %>% filter(species == chosen_capital_species)) %>%
    select(-inversion_ID,-species) %>%
    gather(variable,bayesfactor,-sv_name) 
  associations <- rbind(env_associations,soil_associations)
  if (chosen_species == "petpet" | chosen_species == "petfal"){
    SV_to_keep <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.petinvfreq.txt") %>%
      filter(tolower(species) == chosen_species) %>%
      mutate(chr = gsub("Ha412HOChr","",chr),
             pos=gsub("syn","",mds)) %>%
      select(-mds,-species) %>%
      mutate(sv_name = paste0("pet",chr,".0",pos)) %>%
      filter(freq > min_freq,freq < (1-min_freq)) %>%
      pull(sv_name)
    associations <- associations %>%
      filter(sv_name %in% SV_to_keep)
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
  associations <- associations %>%
    mutate(sv_species_name = paste0(chosen_species,"_",sv_name))
  
  all_associations <- rbind(all_associations,associations)
}
  

variable_type <- read_tsv("gea/environment_variables.txt")

vertical_lines <- c(all_associations %>%distinct(sv_species_name,sv_name) %>%
                      filter(grepl("annuus",sv_species_name)) %>% nrow(),
                    all_associations %>%distinct(sv_species_name,sv_name) %>%
                      filter(grepl("annuus|argophyllus",sv_species_name)) %>% nrow(),
                    all_associations %>%distinct(sv_species_name,sv_name) %>%
                      filter(grepl("annuus|argophyllus|petfal",sv_species_name)) %>% nrow())
bayes_cutoff_1 <- 10
bayes_cutoff_2 <- 20
bayes_cutoff_3 <- 30
header_size <- 5
label_size <- 5
p <- all_associations %>%
  inner_join(.,variable_type) %>%
  mutate(variable = case_when(variable == "sodium" ~ "NA",
                              TRUE ~ variable)) %>%
  ungroup() %>%
  
  mutate(category = case_when(type == "Location" ~ 1,
                              type == "Temperature" ~ 3,
                              type == "Precipitation" ~ 2,
                              type == "Combined" ~ 4,
                              type == "Soil" ~ 5)) %>%
  
  ggplot(.,aes(x=sv_species_name,y=fct_reorder(variable,-category),fill=type)) + 
  geom_tile(alpha=0.2,color="black") +
  scale_fill_brewer(palette="Set1",name="Environmental\ntrait type",
                    breaks=c("Location", "Precipitation", "Temperature","Combined","Soil"),
                    labels=c("Location", "Precipitation", "Temperature","P+T","Soil")) + 
  theme_bw() +
  ylab("Climate variable") +
  xlab("SV") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1,size=label_size),
        axis.text.y = element_text(size=label_size),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8)) +
  scale_x_discrete(labels=all_associations %>% distinct(sv_species_name,sv_name) %>% 
                     arrange(sv_species_name) %>% pull(sv_name)) +
  annotate("text", hjust = 0.5,size=2,
           y=nrow(variable_type),x=(vertical_lines[1]/2)+0.5,label="H. annuus",fontface = 'italic') +
  annotate("text", hjust = 0.5,size=2,
           y=nrow(variable_type),x=(((vertical_lines[2]-vertical_lines[1])/2)+vertical_lines[1])+0.5,label="H. argophyllus",fontface = 'italic') +
  annotate("text", hjust = 0.5,size=2,
           y=nrow(variable_type),x=(((vertical_lines[3]-vertical_lines[2])/2)+vertical_lines[2])+0.5,label="H. p. fallax",fontface = 'italic') +
  annotate("text", hjust = 0.5,size=2,
           y=nrow(variable_type),x=(((all_associations %>% distinct(sv_species_name,sv_name) %>% nrow()-vertical_lines[3])/2)+vertical_lines[3])+0.5,label="H. p. petiolaris",fontface = 'italic') +
  
  coord_cartesian(ylim=c(0.5,nrow(variable_type)))


p <- p +   geom_tile(data=all_associations %>%
                       inner_join(.,variable_type) %>%
                       mutate(variable = case_when(variable == "sodium" ~ "NA",
                                                   TRUE ~ variable)) %>%
                       ungroup() %>%
                       
                       mutate(category = case_when(type == "Location" ~ 2,
                                                   type == "Temperature" ~ 4,
                                                   type == "Precipitation" ~ 3,
                                                   type == "Combined" ~ 1,
                                                   type == "Soil" ~ 5))  %>%
                       filter(bayesfactor > bayes_cutoff_1, bayesfactor <= bayes_cutoff_2),
                     aes(x=sv_species_name,y=fct_reorder(variable,-category)),fill="light grey",alpha=0.8) 
p <- p +   geom_tile(data=all_associations %>%
                       inner_join(.,variable_type) %>%
                       mutate(variable = case_when(variable == "sodium" ~ "NA",
                                                   TRUE ~ variable)) %>%
                       ungroup() %>%
                       
                       mutate(category = case_when(type == "Location" ~ 2,
                                                   type == "Temperature" ~ 4,
                                                   type == "Precipitation" ~ 3,
                                                   type == "Combined" ~ 1,
                                                   type == "Soil" ~ 5))  %>%
                       filter(bayesfactor > bayes_cutoff_2, bayesfactor <= bayes_cutoff_3),
                     aes(x=sv_species_name,y=fct_reorder(variable,-category)),fill="dark grey",alpha=0.8) 
p <- p +   geom_tile(data=all_associations %>%
                       inner_join(.,variable_type) %>%
                       mutate(variable = case_when(variable == "sodium" ~ "NA",
                                                   TRUE ~ variable)) %>%
                       ungroup() %>%
                       
                       mutate(category = case_when(type == "Location" ~ 2,
                                                   type == "Temperature" ~ 4,
                                                   type == "Precipitation" ~ 3,
                                                   type == "Combined" ~ 1,
                                                   type == "Soil" ~ 5))  %>%
                       filter(bayesfactor > bayes_cutoff_3),
                     aes(x=sv_species_name,y=fct_reorder(variable,-category)),fill="black",alpha=0.8) 
p <- p +     geom_vline(xintercept=vertical_lines+0.5,size=2)

print(p)
ggsave("gea/Ha412HO.climategea.sv.baypass.pdf",p, width=183, height=200, units="mm")

