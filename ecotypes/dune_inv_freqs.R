#Plots of all petiolaris inversions with dune ecotype affiliation


library(tidyverse)


sample_info <- read_tsv("sample_info_apr_2018.tsv")
kate_set <- read_tsv("dune_sample_info.txt") %>%
  mutate(location = case_when(spe == "H. neglectus" ~ "MON",
                              spe == "H. petiolaris" ~ "GSD")) %>%
  select(name, location, ecotype)

Nondune_GSD <- c("PET_09")
sample_info %>% filter(population %in% Nondune_GSD) %>%
  select(name) %>% mutate(location = "GSD",ecotype="non-dune") -> nondune_gsd_set

dune_MON <- c("PET_49","PET_47")
sample_info %>% filter(population %in% dune_MON) %>%
  select(name) %>% mutate(location = "MON",ecotype="dune") -> dune_mon_set

nondune_MON <- c("PET_46","PET_48")
sample_info %>% filter(population %in% nondune_MON) %>%
  select(name) %>% mutate(location = "MON",ecotype="non-dune") -> nondune_mon_set

all_dune <- rbind(kate_set, nondune_gsd_set, dune_mon_set,nondune_mon_set)
write_tsv(all_dune, "resources/dune_samples_info.txt")


inv_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt") %>% filter(spe == "petiolaris") %>%
  mutate(LNH = paste0("Ha412HOChr",sprintf("%02d",chr),".",direction,mds))
all_genotypes <- tibble(name=character(), LNH=character(),location=character(),
                        ecotype=character(),state=character(),triangle_genotype=character())
for (i in 1:nrow(inv_list)){
  genotypes <- read_tsv(paste0("MDS_outliers/Ha412HO/petiolaris/Ha412HO_inv.v3.pcasites.Ha412HOChr",
              sprintf("%02d",inv_list$chr[i]),".",inv_list$direction[i],inv_list$mds[i],".genotypes.txt")) %>%
    rename(name = sample) %>%
    mutate(LNH = paste0("Ha412HOChr",sprintf("%02d",inv_list$chr[i]),".",inv_list$direction[i],inv_list$mds[i])) %>%
    inner_join(.,all_dune) %>%
    mutate(state = case_when(location == "GSD" ~"Colorado",
                             location == "MON" ~ "Texas")) %>%
    select(name, LNH,location,ecotype,state,triangle_genotype)
  
  all_genotypes <- rbind(all_genotypes,genotypes)
  
  genotypes <- read_tsv(paste0("MDS_outliers/Ha412HO/petiolaris/Ha412HO_inv.v3.pcasites.Ha412HOChr",
                  sprintf("%02d",inv_list$chr[i]),".",inv_list$direction[i],inv_list$mds[i],".genotypes.txt")) %>%
    rename(name = sample) %>%
    mutate(LNH = paste0("Ha412HOChr",sprintf("%02d",inv_list$chr[i]),".",inv_list$direction[i],inv_list$mds[i])) %>%
    filter(!name %in% all_dune$name) %>%
    inner_join(.,sample_info) %>%
    filter(species == "PetFal" |species =="PetPet") %>%
    mutate(state = "Allopatric",ecotype=species,location="Allopatric") %>%
    select(name, LNH,location,ecotype,state,triangle_genotype) 
  all_genotypes <- rbind(all_genotypes,genotypes)
  
}

all_genotypes %>%
  filter(ecotype != "PetPet" & ecotype != "PetFal") %>%
  group_by(LNH, ecotype, state )%>%
  summarise(gen_freq = mean(triangle_genotype/2)) %>%
  group_by(LNH, state) %>%
  summarize(dif = gen_freq[1]- gen_freq[2]) %>%
  filter(abs(dif) > 0.4) %>%
  select(LNH) -> differential_LNH

pdf("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.duneinversions.pdf",height=4,width=10)
all_genotypes %>%
  inner_join(inv_list) %>%
  filter(ecotype != "PetPet" & ecotype != "PetFal") %>%
  group_by(LNH,sv_name, ecotype, state )%>%
  summarise(gen_freq = mean(triangle_genotype/2)) %>%
  inner_join(differential_LNH) %>%
  separate(.,LNH,c('spacer',"mds"),"HO") %>%
  ggplot(.,aes(x=fct_rev(ecotype),y=gen_freq,color=state)) + 
  geom_hline(data=all_genotypes %>%
               inner_join(inv_list) %>%
               filter(ecotype == "PetPet" | ecotype == "PetFal")%>%
               group_by(LNH,sv_name, ecotype, state )%>%
               summarise(gen_freq = mean(triangle_genotype/2)) %>%
               inner_join(differential_LNH) %>%
               separate(.,LNH,c('spacer',"mds"),"HO"),
             aes(yintercept=gen_freq,linetype=ecotype)) +
  geom_point() +
  geom_line(aes(group=state,color=state)) +
  facet_grid(~sv_name) +
  theme_bw() +
  xlab("") + ylab("Haplotype frequency") +
  scale_color_manual(values=c("#B7A083","#5C4C39"),name="State") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_linetype(name="Subspecies",
                 labels=c("H. petiolaris\nfallax",
                          "H. petiolars\npetiolaris"))

dev.off()




