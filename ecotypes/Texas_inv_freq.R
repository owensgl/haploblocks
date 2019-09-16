#Plots of Texanus allele frequencies

library(tidyverse)
library(SNPRelate)
library(ggrepel)

sample_info <- read_tsv("resources/sample_info_file_all_samples_wildspecies.tsv") %>%
  rename(sample=name,pop=population)
pop_info <- read_tsv("pop_loc_allnum.txt")
texanus_populations <- read_tsv("resources/annuus.texanusfst.grouplist.txt",col_names=c("pop","texanus"))


#PCA of inversions
vcf_out <- "MDS_outliers/Ha412HO/annuus/Ha412HO_inv.v3.pcasites.vcf"
gds_out <- "MDS_outliers/Ha412HO/annuus/Ha412HO_inv.v3.pcasites.gds"
snpgdsVCF2GDS(vcf_out,gds_out, method="biallelic.only",ignore.chr.prefix="Ha412HOChr")
genofile <- snpgdsOpen(gds_out)
pca <- snpgdsPCA(genofile, num.thread=2)
tab <- data.frame(sample = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
snpgdsClose(genofile)

pdf("annuus.texanus.fstplots.pdf",height=5,width=5)
tab %>%
  inner_join(sample_info) %>%
  inner_join(texanus_populations) %>%
  ggplot(.,aes(x=EV1,y=EV2)) + geom_point(aes(color=as.factor(texanus)),alpha=0.9) +
  theme_minimal() +
  scale_color_manual(values=c("#A37212", "#FFD583"),
                     name="Region",labels=c("Central","Texas")) +
  theme(legend.position = "bottom") +
  xlab("PC1") + ylab("PC2")
tab %>%
  inner_join(sample_info) %>%
  inner_join(texanus_populations) %>%
  ggplot(.,aes(EV1)) + geom_density(aes(fill=as.factor(texanus)),alpha=0.5) +
  theme_minimal() +
  scale_fill_manual(values=c("#A37212", "#FFD583"),
                     name="Region",labels=c("Central","Texas")) +
  theme(legend.position = "none") +
  xlab("PC1") +ylab("Density")




#Genome wide PCA, inversions removed
genome_tab <- read_tsv("PCA/Annuus.tranche90.snp.env.90.bi.remappedHa412HO.v3noinv.ldfilter.ldr0p2.pca.txt")
genome_tab %>% 
  rename(sample=name) %>%
  inner_join(sample_info) %>%
  inner_join(texanus_populations) %>%
  ggplot(.,aes(PC1)) + geom_density(aes(fill=as.factor(texanus)),alpha=0.5) +
  theme_minimal() +
  scale_fill_manual(values=c("#A37212", "#FFD583"),
                    name="Region",labels=c("Central","Texas")) +
  theme(legend.position = "none") +
  xlab("PC1") +ylab("Density")

genome_tab %>% 
  rename(sample=name) %>%
  inner_join(sample_info) %>%
  inner_join(texanus_populations) %>%
  ggplot(.,aes(x=PC1,y=PC2)) + geom_point(aes(color=as.factor(texanus)),alpha=0.9) +
  theme_minimal() +
  scale_color_manual(values=c("#A37212", "#FFD583"),
                     name="Region",labels=c("Central","Texas")) +
  theme(legend.position = "bottom")


genome_wide_fst <- read_tsv("/media/owens/Copper/wild_gwas/annuus/Annuus.tranche90.snp.env.90.bi.remappedHa412HO.v3noinv.ldfilter.TEXvCEN.fst.txt.gz")
inv_fst <- read_tsv("/media/owens/Copper/wild_gwas/annuus/Ha412HO_inv.v3.pcasites.TEXvCEN.fst.txt.gz")
inv_fst$percentile <- "NA"
for (i in 1:nrow(inv_fst)){
  n_fst <- inv_fst$Fst[i]
  
  percentile <- genome_wide_fst %>%
    filter(N1 > 5 & N2 > 5) %>%
    mutate(place = case_when(Fst >= n_fst ~ "Higher",
                             TRUE ~ "Lower")) %>%
    group_by(place) %>%
    count() %>%
    ungroup() %>%
    mutate(freq = n / sum(n)) %>%
    filter(place == "Higher") %>%
    pull(freq)
  inv_fst$percentile[i] <- as.numeric(percentile)
}

inv_fst$percent_label <- round(as.numeric(inv_fst$percentile)*100,1)


genome_wide_fst %>%
  ggplot(.,aes(y=Fst,x="Texas")) + 
  geom_violin(fill="#FFD583") +
  geom_point(data=inv_fst,aes(x="Texas",y=Fst)) +
  geom_text_repel(data=inv_fst %>% filter(percent_label < 5),aes(x="Texas",y=Fst,label=paste0(gsub("Ha412HO","", chr),"\n",percent_label,"%")),
                  nudge_x      = 0.2,
                  direction    = "y",
                  angle        = 0,
                  vjust        = 0,
                  segment.size = 0.2
  ) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  xlab("") +
  ylab(expression(F[ST]))

dev.off()

t.test(genome_wide_fst$Fst,inv_fst$Fst)

inversions <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt") %>%
  filter(spe == "annuus")


all_genos <- tibble(sample=character(),triangle_genotype=numeric(),pop=character(),species=character(),
                    inversion=character())
for (i in 1:11){
  

chr_data <- read_tsv(paste0("MDS_outliers/Ha412HO/annuus/Ha412HO_inv.v3.pcasites.Ha412HOChr",
                            sprintf("%02d",inversions$chr[i]),".",
                            inversions$direction[i],inversions$mds[i],".genotypes.txt"))

usa <- map_data('state')
states <- map_data("state")
target_state <- map_data('state')
lat_range <- c(25, 50)
long_range <- c(-125,-93)
pie_size <- 0.4


chr_data %>%
  inner_join(.,sample_info) %>%
  filter(species== "Ann") %>%
  inner_join(.,pop_info) %>%
  group_by(pop,lat,long) %>%
  summarize(mean_loci=mean(triangle_genotype)) -> plotting_data

chr_data %>%
  inner_join(.,sample_info) %>%
  filter(species== "Ann") %>%
  inner_join(.,pop_info) %>%
  mutate(inversion = paste0("Chr",sprintf("%02d",inversions$chr[i]),".",
                            inversions$direction[i],inversions$mds[i])) %>%
  select(sample,triangle_genotype,pop,species,inversion) -> tmp
all_genos <- rbind(all_genos,tmp)
  

}
all_genos %>%
  inner_join(pop_info) %>%
  group_by(inversion,pop,lat,long) %>%
  summarize(mean_loci=mean(triangle_genotype)) -> plotting_data

ggplot(target_state, aes(long, lat)) +
  geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
  coord_quickmap() +
  geom_point(data=plotting_data,
             aes(x=long, y=lat, color=mean_loci/2), 
             size=3,alpha=0.8) +
  scale_color_viridis_c(name="Allele_frequency") +theme_bw() +
  xlab("Longitude") + ylab("Latitude") +
  scale_x_continuous(limits = long_range, expand = c(0, 0)) +
  scale_y_continuous(limits = lat_range, expand = c(0, 0)) +
  facet_wrap(~inversion)

all_genos %>%
  inner_join(.,pop_info) %>%
  inner_join(.,texanus_populations) %>%
  group_by(inversion,texanus) %>%
  summarize(mean_loci=mean(triangle_genotype)) %>%
  ggplot(.,aes(x=texanus,y=mean_loci/2)) + geom_bar(stat="identity",position="dodge") +
  facet_wrap(~inversion)
