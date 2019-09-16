#This is for making the PCA plot of the argophyllus chr6 FT region
library(tidyverse)
library(SNPRelate)
het_script <- "/home/owens/bin/pop_gen/vcf2het.pl"

data_directory <- "/media/owens/Copper/wild_gwas/"
species_directory <- "argophyllus"
species <- "Argophyllus"
vcf_tag <- "gwas"
vcf_in <- paste(
  data_directory,
  "/",
  species_directory,
  "/",
  species,
  ".tranche90.snp.",
  vcf_tag,
  ".90.bi.remappedHa412HO.vcf.gz",
  sep="")
vcf_out <- paste(
  data_directory,
  "/",
  species_directory,
  "/",
  species,
  ".tranche90.snp.",
  vcf_tag,
  ".90.bi.remappedHa412HO.tmp.vcf.gz",
  sep="")
gds_out <- paste(
  data_directory,
  "/",
  species_directory,
  "/",
  species,
  ".tranche90.snp.",
  vcf_tag,
  ".90.bi.remappedHa412HO.tmp.gds",
  sep="")
region_string <- "Ha412HOChr06:130156472-155128857"


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

try_3_clusters <-try(kmeans(tab[,2], 3, centers=c(min(tab[,2]),(min(tab[,2])+max(tab[,2]))/2,max(tab[,2]))))


if("try-error" %in% class(try_3_clusters)){
  kmeans_cluster <-kmeans(tab[,2], 2, centers=c(min(tab[,2]),max(tab[,2])))
  
}else{
  
  kmeans_cluster <-kmeans(tab[,2], 3, centers=c(min(tab[,2]),(min(tab[,2])+max(tab[,2]))/2,max(tab[,2])))
  
}

tab$cluster <- kmeans_cluster$cluster - 1 

betweenss <- kmeans_cluster$betweenss

max_het <- round(max(vcf_het$percent_het),2)
min_het <- round(min(vcf_het$percent_het),2)
pdf("figures/ID19_Arg_chr6PCA.pdf")
tab %>% 
  ggplot(.,aes(x=EV1,y=EV2)) + geom_point(aes(color=as.factor(cluster)),size=2) +
  theme_bw() + 
  scale_color_manual(values=paper_colors,name="Kmeans\ncluster") +
  scale_shape(name="Kmeans\ncluster") +
  theme(legend.position="bottom") +
  ylab("PC2") + xlab("PC1")
ggsave("figures/ID19_Arg_chr6PCA.pdf",device="pdf")

tab %>%
  inner_join(vcf_het) %>%
  ggplot(.,aes(x=as.factor(cluster),y=percent_het)) + geom_boxplot()

paper_colors <- c("light grey","dark grey","black")
vcf_het %>%
  inner_join(tab) %>%
  ggplot(.,aes(x=as.factor(cluster),y=percent_het)) + 
  geom_violin(aes(fill=as.factor(cluster),color=as.factor(cluster))) +
  scale_fill_manual(values=paper_colors,name="Kmeans\ncluster") +
  scale_color_manual(values=paper_colors,name="Kmeans\ncluster") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Heterozygosity") +
  labs(tag = "C") +
  theme(legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
ggsave("figures/ID68_Arg_chr6PCA_heterozygosity.pdf",device="pdf")

