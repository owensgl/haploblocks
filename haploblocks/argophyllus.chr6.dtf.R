#Correlation between FT and inversion genotypes across ARG chr6 inversions. 

library(tidyverse)
library(ggthemes)
library(forcats)
library(ggrastr)
library(gridExtra)
library(broom)
library(tidyquant)
library(zoo)
library(beeswarm)



#Load in all inversion genotypes
pca_genotypes <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.inversions.genotypes.v1.txt")
folder <- "MDS_outliers"
chosen_species <- "argophyllus"
chosen_species_abbreviation <- "Arg"
chosen_abbreviation <- "ARG"
filtering <- "pcasites"
prefix <- "Ha412HO_inv.jan09"
base_directory <- "/media/owens/Copper/wild_gwas/"
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T)
inversions <- pca_genotypes %>% filter(species == chosen_species) %>%
  select(chr,mds) %>% filter(chr == "Ha412HOChr06") %>%unique()

#Load the genotype data for the two chr6 inversions.
full_genotypes <- tibble(chr=character(),pos=numeric(),sample=character(),genotype=numeric(),
                         depth=numeric())
for (n in 1:2){
  chosen_chr <- inversions[n,1] %>% pull()
  chosen_mds <- inversions[n,2] %>% pull()
  min_depth_genotyping <- 2
  
  fst_calls <- read_tsv(paste(base_directory,chosen_species,"/",prefix,".",chosen_species,".", chosen_chr,".",chosen_mds,".",
                              filtering,".genotypes.txt.gz",sep=""),guess_max= 100000) %>%
    inner_join(., labels %>% dplyr::rename(sample=name)) %>%
    filter(genotype == "00" | genotype == "01" | genotype == "11") %>%
    mutate(genotype = case_when(genotype == "00" ~ 0,
                                genotype == "01" ~ 1,
                                genotype == "11" ~ 2)) %>%
    filter(species == chosen_species_abbreviation) %>%
    filter(depth >= min_depth_genotyping)
  full_genotypes <- rbind(full_genotypes,fst_calls %>% select(chr,pos,sample,genotype,depth))
}

phenotypes <- read_tsv(paste("resources/ARG_all_phenotypes_May2019.txt",sep="")) %>%
  rename(sample = FID) %>% select(sample, DTF_deduced) %>%
  mutate(early_late = case_when(DTF_deduced < 80 ~ "Early",
                                TRUE ~ "Late"))


full_genotypes <- full_genotypes %>% inner_join(phenotypes)





  
#pdf("Ha412HO_inv.jan09.inversions.argophyllus.late_early.pdf",height=20,width=8)



full_genotypes %>%
  select(pos) %>% unique() %>% pull() -> positions

n_labels <- ceiling(length(positions)/20)
x_labels <- sort(positions)[seq(0, length(positions), by= n_labels)]
full_genotypes %>%
  group_by(sample) %>%
  mutate(sum_geno =sum(genotype)) %>%
  ungroup() %>%
  mutate(sample = as.factor(sample)) %>%
  mutate(sample = fct_reorder(sample, early_late)) %>% 
  #complete(sample,nesting(pos)) %>%
  group_by(sample) %>%
  arrange(sample, pos) %>%
 # mutate(mean5genotype = rollmean(x = genotype, 10, align = "right", fill = NA,na.rm=T)) %>%
  ggplot(.,aes(x=as.factor(pos),y=fct_reorder(sample, sum_geno))) + geom_tile(aes(fill=as.factor(genotype))) +
    scale_fill_manual(name="Genotype",values=c("light grey","dark grey","black")) +
    facet_wrap(~early_late,nrow=2,scales="free_y") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.text.x=element_text(angle=60, hjust=1)) +
  scale_x_discrete(breaks=x_labels)  +
  xlab("Position") +
  ylab("Sample")
  
ggsave("Ha412HO_inv.jan09.inversions.argophyllus.late_early.png",units="in",width=8,height=4,device="png")
  
dev.off()

#Plotting the flowering time

genotypes <- read_tsv("MDS_outliers/Ha412HO/argophyllus/Ha412HO_inv.v3.pcasites.Ha412HOChr06.10MB.genotypes.txt",
                      col_names = c("sample","genotype"))
arg_ft <- read_tsv("gwas/argophyllus/DTF.txt",col_names = c("sample","id","dtf"))

genotypes %>%
  inner_join(arg_ft) %>%
  mutate(genotype = case_when(genotype == "1" ~ 2,
                              genotype == "0" ~ 0,
                              genotype == "2" ~ 2)) -> dominant_genotypes
  anova(lm(dtf ~ genotype ,dominant_genotypes))
  af <- anova(lm(dtf ~ genotype ,dominant_genotypes))
  afss <- af$"Sum Sq"
  print(cbind(af,PctExp=afss/sum(afss)*100))

ann_ft <- read_tsv("gwas/annuus/DTF.txt",col_names = c("sample","id","dtf")) %>%
  mutate(genotype = "ANN")

pdf("Ha412HO_inv.v3.pcasites.Ha412HOChr06.10MB.genotypes.dtf.pdf",height=4,width=4)
arg_ft %>%
  inner_join(genotypes) %>%
  rbind(.,ann_ft) %>%
  ggplot(.,aes(x=genotype,y=dtf,color=genotype)) + 
  #geom_boxplot(color="black") +
  geom_jitter(width=0.2,alpha=0.7) +
  scale_color_manual(values=c("#6DBF76","#3DA047","#0D6D16","#FFB31C")) +
  theme_minimal() +
  ylab("Days to flower") +
  theme(legend.position = "none")

arg_ft %>%
  inner_join(genotypes) %>%
  rbind(.,ann_ft) %>%
  ggplot(.,aes(x=genotype,y=dtf,color=genotype)) + 
  #geom_boxplot(color="black") +
  ggbeeswarm::geom_quasirandom(width=0.4,alpha=0.7) +
  scale_color_manual(values=c("#6DBF76","#3DA047","#0D6D16","#FFB31C")) +
  theme_minimal() +
  ylab("Days to flower") +
  theme(legend.position = "none")

arg_ft %>%
  inner_join(genotypes) %>%
  rbind(.,ann_ft) %>%
  ggplot(.,aes(x=genotype,y=dtf,fill=genotype)) + 
  geom_boxplot() +
  #ggbeeswarm::geom_quasirandom(width=0.2,alpha=0.8) +
  scale_fill_manual(values=c("#6DBF76","#3DA047","#0D6D16","#FFB31C")) +
  theme_minimal() +
  ylab("Days to flower") +
  theme(legend.position = "none")
dev.off()

####Heterozygosity of chr6 inversion samples

het <- read_tsv("MDS_outliers/Ha412HO/argophyllus/Annuus.tranche90.snp.annarg.90.bi.chr6test.remappedHa412HO.het.txt")
sample_types <- read_delim("MDS_outliers/Ha412HO/argophyllus/chr6_invtest.sampleinfo.txt",col_names = c("sample","type"),delim=" ")


pdf("MDS_outliers/Ha412HO/argophyllus/Annuus.tranche90.snp.annarg.90.bi.chr6test.remappedHa412HO.het.pdf")
het %>%
  inner_join(.,sample_types) %>%
  filter(type != "ann") %>%
  mutate(type = case_when(type == "early" ~ "Early flowering",
                          TRUE ~ "Late flowering")) %>%
  ggplot(.,aes(x=type,y=percent_het)) + geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  ylab("Heterozygosity") +
  theme(axis.title.x=element_blank())
dev.off()
  
  
##Looking at depth in inversion area
phenotypes <- read_tsv(paste("/media/owens/Copper/wild_gwas_2018/gwas/",chosen_species,"/",chosen_abbreviation,"_all_phenotypes_jan_2019.txt",sep="")) %>%
  rename(sample = FID) %>% select(sample, DTF) %>%
  mutate(DTF_group = case_when(DTF< 100 ~ "Early",
                               TRUE ~ "Late"))

depths <- read_tsv("/media/owens/Copper/wild_gwas_2018/argophyllus/argophyllus.chr6.tidydepth.txt.gz")
  
depths %>%
  inner_join(.,phenotypes) -> depths

window_size = 100000
depths %>% 
  mutate(window = floor(pos/window_size)*window_size,
         present = case_when(depth == 0 ~ 0,
                             TRUE ~ 1)) %>%
  group_by(window,sample) %>%
  mutate(sites = n()) %>%
  group_by(DTF_group,window,chr) %>%
  summarize(sequenced = mean(present), sites = mean(sites)) -> processed_depth

processed_depth %>%
  filter(sites > 10) %>%
  group_by(chr,window) %>%
  summarize(sequenced_dif = (sequenced[which(DTF_group == "Early")] - 
              sequenced[which(DTF_group == "Late")])) %>%
  ggplot(.,aes(x=window,y=sequenced_dif)) + geom_line()


processed_depth %>%
  filter(sites > 10) %>%
  #filter(DTF_group == "Late") %>%
  ggplot(.,aes(x=window,y=sequenced,color=DTF_group)) + geom_line() +
  theme_bw()


depths %>%
  mutate(present = case_when(depth == 0 ~ "0",
                             TRUE ~ "1")) %>%
  filter(pos > 134214322, pos <134219438) %>% 
  #filter(pos > 132000000, pos <135000000) %>%
  mutate(sample = as.factor(sample)) %>%
  ggplot(.,aes(x=as.factor(pos),y=fct_reorder(sample,DTF),fill=present)) + geom_tile() +
  scale_fill_manual(values=c("#ffffff","#000000"))

         