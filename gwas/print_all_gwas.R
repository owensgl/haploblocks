#Plotting all Inversion GWAS scores 

library(tidyverse)
library(RColorBrewer)
library(ggrastr)

pvalue_cutoff <- 0.001
min_freq <- 0.03
directory <- "/media/owens/Copper/wild_gwas/gwas"
species_abbreviation_list <- c("annuus","argophyllus","petpet","petfal")
species_capital_list <-c("Annuus", "Argophyllus","Petiolaris","Petiolaris")
tag_list <- c("gwas", "gwas","petpet","petfal")
trait_types <- read_tsv("gwas/all_traits_category.txt",col_names = c("trait","trait_type"))

inversion_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt") %>%
  mutate(inversion_ID = paste0("Ha412HOChr",sprintf("%02d",chr),"_",direction,mds)) %>%
  select(inversion_ID, sv_name,species) %>%
  mutate(id = inversion_ID) %>%
  mutate(id = gsub("Ha412HOChr","",id)) %>%
  mutate(id = gsub("pos","",id)) %>%
  mutate(id = gsub("neg","",id)) %>%
  mutate(id = gsub("syn","",id))

all_data <- tibble(sv_name=character(),sv_species_name=character(),trait=character(),pvalue=numeric())
#pdf("gwas/Ha412HO_inv.v3.gwas.threshold.pdf",height=14,width=12)
for (n in 1:4){
  
  
  species_abbreviation <- species_abbreviation_list[n]
  species_capital <- species_capital_list[n]
  tag <- tag_list[n]
  
  ##Get list of all traits
  traits <- tibble(trait=list.files(paste0("gwas/",species_abbreviation),include.dirs = F) ) %>%
    filter(!grepl("znorm",trait)) %>%
    separate(trait,c("trait","spacer"),"\\.") %>%
    select(-spacer)
  
  
  #Pull out associations for inversions
  
  all_invs <- read_tsv(paste0(directory,"/",species_abbreviation,"/results/",
                              "Ha412HO_inv.v3.pcasites.",species_abbreviation,
                              ".shuffled.",
                              traits$trait[1],".ps.gz"),
                       col_names=c("id","beta","sd","pvalue")) %>%
    select(id) %>% separate(id,c("chr","pos"),"_") %>%
    mutate(count=0)
  all_gwas_inv <- read_tsv(paste0(directory,"/",species_abbreviation,"/results/",
                                               "Ha412HO_inv.v3.pcasites.",species_abbreviation,
                                               ".shuffled.",
                                               traits$trait[1],".ps.gz"),
                                        col_names=c("id","beta","sd","pvalue")) %>% mutate(trait = traits$trait[1],
                                                                                           species = species_abbreviation) %>%head(0)
  for (x in 1:nrow(traits)){
    gwas_inv <- read_tsv(paste0(directory,"/",species_abbreviation,"/results/",
                                 "Ha412HO_inv.v3.pcasites.",species_abbreviation,
                                 ".shuffled.",
                                 traits$trait[x],".ps.gz"),
                          col_names=c("id","beta","sd","pvalue")) %>%
      mutate(trait = traits$trait[x],
             species = species_abbreviation)
    all_gwas_inv <- rbind(all_gwas_inv,gwas_inv)

    

  }
  all_gwas_inv <- inner_join(all_gwas_inv, inversion_list %>% filter(species == species_capital) %>% select(id,sv_name)) %>%
    mutate(sv_species_name = paste0(species,"_",sv_name))
  
  ##Remove inversions below a minimum allele frequency
  if (species_capital == "Petiolaris"){
    inv_frequencies <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.petinvfreq.txt") %>%
      filter(tolower(species) == species_abbreviation) %>%
      mutate(chr = gsub("Ha412HOChr","",chr),
             pos=gsub("syn","",mds)) %>%
      select(-mds,-species) %>%
      mutate(id = paste0(chr,"_",pos))
    inner_join(all_gwas_inv, inv_frequencies) %>%
      filter(freq > min_freq,freq < (1-min_freq)) %>%
      select(-freq) -> all_gwas_inv

  }
  all_gwas_inv <- all_gwas_inv %>%
    select(sv_name,sv_species_name,trait,pvalue) 
  
  all_data <- rbind(all_data,all_gwas_inv)
}
#Make figure

trait_types <- trait_types %>% filter(trait %in% unique(all_data$trait))
vertical_lines <- c(all_data %>%distinct(sv_species_name,sv_name) %>%
                      filter(grepl("annuus",sv_species_name)) %>% nrow(),
                    all_data %>%distinct(sv_species_name,sv_name) %>%
                      filter(grepl("annuus|argophyllus",sv_species_name)) %>% nrow(),
                    all_data %>%distinct(sv_species_name,sv_name) %>%
                      filter(grepl("annuus|argophyllus|petfal",sv_species_name)) %>% nrow())
pvalue_cutoff_1 <- 0.01
pvalue_cutoff_2 <- 0.001
pvalue_cutoff_3 <- 0.0001
label_size <- 5
p <- all_data %>%
  inner_join(.,trait_types) %>%
  ungroup() %>%
  mutate(significant = case_when(pvalue < pvalue_cutoff ~ "+",
                                 TRUE ~ "")) %>% 
  mutate(category = case_when(trait_type == "Architecture" ~ 1,
                              trait_type == "CN_content" ~ 2,
                              trait_type == "Development" ~ 3,
                              trait_type == "Inflorescence" ~ 4,
                              trait_type == "Leaf" ~ 5,
                              trait_type == "Pigmentation" ~ 6,
                              trait_type == "Seed" ~ 7)) %>%
  
  ggplot(.,aes(x=sv_species_name,y=fct_reorder(trait,-category),fill=trait_type)) + 
  geom_tile(alpha=0.2,color="black") +
  scale_fill_manual(values=brewer.pal(7,"Set1"),name="Trait type") + 
  #geom_tile(aes(label=significant)) +
  theme_bw() +
  ylab("Trait") +
  xlab("SV") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1,size=label_size),
        axis.text.y = element_text(size=label_size),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        legend.text=element_text(size=6)) +
  scale_x_discrete(labels=all_data %>% distinct(sv_species_name,sv_name) %>% 
                     arrange(sv_species_name) %>% pull(sv_name)) +
  annotate("text", hjust = 0.5,size=2,
           y=nrow(trait_types)+1,x=(vertical_lines[1]/2)+0.5,label="H. annuus",fontface = 'italic') +
  annotate("text", hjust = 0.5,size=2,
           y=nrow(trait_types)+1,x=(((vertical_lines[2]-vertical_lines[1])/2)+vertical_lines[1])+0.5,label="H. argophyll.",fontface = 'italic') +
  annotate("text", hjust = 0.5,size=2,
           y=nrow(trait_types)+1,x=(((vertical_lines[3]-vertical_lines[2])/2)+vertical_lines[2])+0.5,label="H. p. fallax",fontface = 'italic') +
  annotate("text", hjust = 0.5,size=2,
           y=nrow(trait_types)+1,x=(((all_data %>% distinct(sv_species_name,sv_name) %>% nrow()-vertical_lines[3])/2)+vertical_lines[3])+0.5,label="H. p. petiolaris",fontface = 'italic') +
  
  coord_cartesian(ylim=c(0.5,nrow(trait_types)+1))

p <- p +   geom_tile(data=all_data %>%
                       inner_join(.,trait_types) %>%
                       ungroup() %>%
                       mutate(category = case_when(trait_type == "Architecture" ~ 1,
                                                   trait_type == "CN_content" ~ 2,
                                                   trait_type == "Development" ~ 3,
                                                   trait_type == "Inflorescence" ~ 4,
                                                   trait_type == "Leaf" ~ 5,
                                                   trait_type == "Seed" ~ 6)) %>%
                       filter(pvalue < pvalue_cutoff_1, pvalue >= pvalue_cutoff_2),
                     aes(x=sv_species_name,y=fct_reorder(trait,-category)),fill="light grey",alpha=0.8) 
p <- p +   geom_tile(data=all_data %>%
                       inner_join(.,trait_types) %>%
                       ungroup() %>%
                       mutate(category = case_when(trait_type == "Architecture" ~ 1,
                                                   trait_type == "CN_content" ~ 2,
                                                   trait_type == "Development" ~ 3,
                                                   trait_type == "Inflorescence" ~ 4,
                                                   trait_type == "Leaf" ~ 5,
                                                   trait_type == "Seed" ~ 6)) %>%
                       filter(pvalue < pvalue_cutoff_2, pvalue >= pvalue_cutoff_3),
                     aes(x=sv_species_name,y=fct_reorder(trait,-category)),fill="dark grey",alpha=0.8) 
p <- p +   geom_tile(data=all_data %>%
                       inner_join(.,trait_types) %>%
                       ungroup() %>%
                       mutate(category = case_when(trait_type == "Architecture" ~ 1,
                                                   trait_type == "CN_content" ~ 2,
                                                   trait_type == "Development" ~ 3,
                                                   trait_type == "Inflorescence" ~ 4,
                                                   trait_type == "Leaf" ~ 5,
                                                   trait_type == "Seed" ~ 6)) %>%
                       filter(pvalue < pvalue_cutoff_3),
                     aes(x=sv_species_name,y=fct_reorder(trait,-category)),fill="black",alpha=0.8) 
p <- p +     geom_vline(xintercept=vertical_lines+0.5,size=2)

print(p)

ggsave("gwas/Ha412HO_inv.v3.gwas.threshold.pdf",p, width=183, height=200, units="mm")




##########Print all GWAS.

chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")
for (n in 4:4){
  
  
  species_abbreviation <- species_abbreviation_list[n]
  species_capital <- species_capital_list[n]
  tag <- tag_list[n]
  pdf(paste0("gwas/Ha412HO_inv.v3.gwas.",species_abbreviation,".withoutSV.SNPs.pdf"),height=6,width=12)
  
  if (species_abbreviation == "annuus"){
    species_full_name <- "H. annuus"
  }else if (species_abbreviation == "argophyllus"){
    species_full_name <- "H. argophyllus"
  }else if (species_abbreviation == "petpet"){
    species_full_name <- "H. petiolaris petiolaris"
  }else if (species_abbreviation == "petfal"){
    species_full_name <- "H. petiolaris fallax"
  }
  
  ##Get list of all traits
  traits <- tibble(trait=list.files(paste0("gwas/",species_abbreviation),include.dirs = F) ) %>%
    filter(!grepl("znorm",trait)) %>%
    separate(trait,c("trait","spacer"),"\\.") %>%
    select(-spacer)
  freqs <- read_tsv(paste0("/media/owens/Copper/wild_gwas/gwas/",species_abbreviation,"/",species_capital,".tranche90.snp.",tag,".90.bi.remappedHa412HO.beagle.freq"),
                    col_names=c("chr","pos","freq")) %>%
    filter(freq > min_freq, freq < (1-min_freq)) %>%
    select(-freq)
  

  for (x in 1:nrow(traits)){
    chosen_trait <- traits$trait[x]
    
    
    inv_associations <- read_tsv(paste0("/media/owens/Copper/wild_gwas/gwas/",species_abbreviation,"/results/","Ha412HO_inv.v3.pcasites.",species_abbreviation,
                                                 ".shuffled.",chosen_trait,".ps.gz"),
                                 col_names=c("id","beta","sd","pvalue")) %>%
      mutate(variable = chosen_trait,
             logp = abs(log10(pvalue)))
    min_freq <- 0.03
    if (species_abbreviation == "petpet" | species_abbreviation == "petfal"){
      inv_frequencies <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.petinvfreq.txt") %>%
        filter(tolower(species) == species_abbreviation) %>%
        mutate(chr = gsub("Ha412HOChr","",chr),
               pos=gsub("syn","",mds)) %>%
        select(-mds,-species) %>%
        mutate(id = paste0(chr,"_",pos))
      inner_join(inv_associations, inv_frequencies) %>%
        filter(freq > min_freq,freq < (1-min_freq)) %>%
        select(-freq) -> inv_associations
    }

    #Load up inversion locations
    inv_locations <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt") %>%
      mutate(full_id = paste0(chr,"_",mds)) %>%
      mutate(id = gsub("Ha412HOChr","",full_id)) %>%
      mutate(id = gsub("pos","",id)) %>%
      mutate(id = gsub("neg","",id)) %>%
      mutate(id = gsub("syn","",id)) %>%
      filter(species == species_capital)
    
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
    

    
    gwas <- read_tsv(paste0("/media/owens/Copper/wild_gwas/gwas/",species_abbreviation,"/results/",species_capital,".tranche90.snp.",tag,
                            ".90.bi.remappedHa412HO.",chosen_trait,".beagle.fullgenome.v3noinv.ldfilter.ps.gz"),
                     col_names=c("id","beta","sd","pvalue")) %>%
      separate(id,c("chr","pos"),"_") %>%
      mutate(pos=as.numeric(pos)) %>%
      semi_join(.,freqs) %>%
      mutate(chr_n = chr, chr= paste0("Ha412HOChr",chr_n))


    gwas_cum <- chr_lengths %>%
      select(chr,end) %>%
      rename(chr_len = end) %>%
      
      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(chr_len)-chr_len) %>%
      select(-chr_len) %>%
      
      # Add this info to the initial dataset
      left_join(gwas, ., by=c("chr"="chr")) %>%
      
      # Add a cumulative position of each SNP
      arrange(chr, pos) %>%
      mutate( poscum=pos+tot)


    axisdf = gwas_cum %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )

    inversion_to_be_plotted <- inv_associations %>%
      inner_join(.,inv_locations)
    
    p <- gwas_cum %>%
      mutate(logp = abs(log10(pvalue))) %>%
      filter(logp > 2) %>%
      ggplot(.) +
      geom_point_rast( aes(x=poscum, y=logp,color=as.factor(chr)), alpha=0.8, size=1.3) +
      # Show all points
      scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +

      # custom X axis:
      scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
      theme_bw() +
      theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
      ggtitle(paste(species_full_name,"-",chosen_trait,"without SV regions")) +
      xlab("Chr") +
      geom_segment(data=inversion_to_be_plotted %>%
                     filter(as.numeric(logp) > 2),
                   aes(x=startcum,xend=endcum,
                       y=as.numeric(logp),yend=as.numeric(logp)),
                   size=2,color="black",alpha=0.8)
    
    print(p)
  }
  dev.off()
  
}

# pdf(paste0("gwas/Ha412HO_inv.v3.gwas.",species_abbreviation,".exampleplateau.SNPs.pdf"),height=4,width=12)
# 
# gwas_cum %>%
#   mutate(logp = abs(log10(pvalue))) %>%
#   filter(logp > 2) %>%
#   ggplot(.) +
#   geom_point_rast( aes(x=poscum, y=logp,color=as.factor(chr)), alpha=0.8, size=1.3) +
#   # Show all points
#   scale_color_manual(values = rep(c("grey", "#447499"), 22 )) +
#   
#   # custom X axis:
#   scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
#   theme_bw() +
#   theme(legend.position="none",
#         panel.border = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank()) +
#   xlab("Chr") + ylab("-log10 p-value") +
#   geom_hline(yintercept=8,linetype="dashed")
# dev.off()


# 
# axisdf = gwas_cum %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )
# 
# 
# 
# 
# sample_info <- read_tsv("resources/sample_info_file_all_samples_2018_12_07.tsv") %>% rename(pop = population,sample=name)
# pop_loc <- read_tsv("pop_loc_allnum.txt")
# genotypes <- read_tsv("MDS_outliers/Ha412HO/annuus/Ha412HO_inv.v3.pcasites.Ha412HOChr13.pos1.genotypes.txt")
# traits <- read_tsv("gwas/annuus/Days_to_budding.txt",col_names = c("sample","spacer","DTB"))
# pca <- read_delim("PCA/Annuus.tranche90.snp.gwas.90.bi.remappedHa412HO.ldr0p2.eigenvectors.txt",delim=" ",
#                   col_names = paste0("PCA",seq(1:614)))
# cbind(traits %>% select(sample),pca ) %>% as_tibble() -> pca
# af <- anova(lm(DTB ~ PCA1 + PCA2 +PCA3+triangle_genotype, 
#                genotypes %>%
#                  inner_join(traits) %>%
#                  inner_join(pca) %>%
#                  inner_join(sample_info)
#                  #filter(pop != "ANN_48", pop != "ANN_49")
#                  ))
# 
# afss <- af$"Sum Sq"
# print(cbind(af,PctExp=afss/sum(afss)*100))
# 
# pca %>%
#   inner_join(sample_info) %>% inner_join(pop_loc) %>%
#   inner_join(traits) %>% 
#   inner_join(genotypes) %>%
#   group_by(pop,lat,long) %>% 
#   group_by(triangle_genotype) %>%
#   filter(pop != "ANN_48", pop != "ANN_49") %>%
#   summarize(mean= mean(DTB,na.rm=T))
#   ggplot(.,aes(x=as.factor(triangle_genotype),y=DTB)) + geom_jitter(aes(color=DTB),size=2) +scale_color_viridis() + theme_minimal()
# 
# 
# #Try individual sites
# samples <- read_tsv("/media/owens/Copper/wild_gwas/annuus/Annuus.tranche90.snp.gwas.90.bi.samplelist.txt",col_names=c("sample")) 
# system(paste("bcftools query -H -f '%END [ %GT]\n' -r Ha412HOChr13:9045002 /media/owens/Copper/wild_gwas/annuus/Annuus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.bcf",
#              ' | sed s/\\#\\ //g |  sed s/\\:GT//g | sed s/END/pos/g > tmp.geno.txt',sep=""))
# 
# read_delim("tmp.geno.txt",delim=" ",col_names = c("pos",as.vector(samples$sample)),skip=1) %>%
#   mutate_if(.,
#                                is.character,
#                                str_replace_all, pattern = "0\\|0", replacement = "0") %>%
#   mutate_if(.,
#             is.character,
#             str_replace_all, pattern = "1\\|1", replacement = "2") %>%
#   mutate_if(.,
#             is.character,
#             str_replace_all, pattern = "0\\|1", replacement = "1") %>%
#   mutate_if(.,
#             is.character,
#             str_replace_all, pattern = "1\\|0", replacement = "1") %>%
#   gather(sample,genotype,ANN0801:ANN1501)%>%
#   mutate(genotype = as.numeric(genotype))-> snps
#          
# snps %>%
#   inner_join(.,traits) %>%
#   inner_join(pca) 
#   
# af <- anova(lm(DTB ~ PCA1 + PCA2 +PCA3+genotype, 
#                snps %>%
#                  inner_join(.,traits) %>%
#                  inner_join(pca) ))
# 
# afss <- af$"Sum Sq"
# print(cbind(af,PctExp=afss/sum(afss)*100))
# q  
#   
  
