library(tidyverse)

#Plotting missing data in ARG chr6 region.

pure_early <- c("ARG0180","ARG0177","ARG0143","ARG0132","ARG0122","ARG0172","ARG0171","ARG0140","ARG0128","ARG0166","ARG0188","ARG0161")
chr6_genotypes <- read_tsv("MDS_outliers/Ha412HO/argophyllus/Ha412HO_inv.jan09.pcasites.Ha412HOChr06.pos1.genotypes.txt")

chr6_genotypes %>%
  filter(triangle_genotype == "0",perc_1  < 0.05) %>% pull(sample) -> pure_late

reads <- read_tsv("/media/owens/Copper/wild_gwas/argophyllus/argophyllus.chr6.tidydepth.txt.gz")
reads %>% filter(sample %in% pure_early | sample %in% pure_late) %>% 
  filter(pos > 130000000, pos < 140000000) %>%
  mutate(group = case_when(sample %in% pure_early ~ "early",
                            sample %in% pure_late ~ "late")) -> filtered_reads

window_size = 50000
filtered_reads  %>%
  mutate(any_depth = case_when(depth == 0 ~ 0,
                               depth > 0 ~ 1)) %>%
  group_by(pos,group) %>%
  summarize(perc_present = mean(any_depth,na.rm=T)) %>%
  mutate(window = (floor(pos/window_size)*window_size)+(window_size/2),
         site_missing = case_when(perc_present < 0.05 ~ 0,
                   TRUE ~ 1)) %>%
  group_by(window,group) %>%
  summarize(percent_missing_sites = mean(site_missing),
            n_sites=n()) -> plotting_data
  
plotting_data %>%
  ggplot(.,aes(x=window,y=percent_missing_sites,color=percent_missing_sites)) + geom_point() +
  facet_wrap(~group,nrow=2) + theme_bw()
