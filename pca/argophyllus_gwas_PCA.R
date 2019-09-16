library(SNPRelate)
library(tidyverse)
library(viridis)
library(gridExtra)
library(purrr)

folder <- "/media/owens/Copper/wild_gwas/argophyllus/"
gds.file <- "Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.gds"
vcf.file <- "Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.vcf.gz"
snpgdsVCF2GDS(paste(folder,vcf.file,sep="/"), paste(folder,gds.file,sep="/"), method="biallelic.only", ignore.chr.prefix = "HanXRQChr")

genofile <- snpgdsOpen(paste(folder,gds.file,sep="/"))
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,method="r", num.thread=10,autosome.only = F)
snpset.id <- unlist(snpset)

pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=10, 
                 eigen.cnt = 0,autosome.only = F)

#Write eigenvectors and eigenvalues
write(pca$eigenval, "PCA/Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.ldr0p2.eigenvalues.txt",
      ncol=length(pca$eigenval))


write.table(pca$eigenvect, "PCA/Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.ldr0p2.eigenvectors.txt",
            col.names = F,row.names = F)


pc.percent <- pca$varprop*100
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T)
pop_loc <- read_tsv("pop_loc_allnum.txt")
pop_loc %>% rename(population = pop) %>% inner_join(.,labels) -> labels

snpgdsClose(genofile)


#Plotting PCAs
pdf("Argophyllus.tranche90.snp.gwas.90.bi.ldr0p2.pca.pdf",height=6,width=14)
for (i in seq(1,20,2)){
  j <- i +1
  tab <- data.frame(name = pca$sample.id,
                    EV1 = pca$eigenvect[,i],  
                    EV2 = pca$eigenvect[,j],    
                    
                    stringsAsFactors = FALSE)
  inner_join(tab, labels) -> tab
  lat <- ggplot(tab) + geom_point(aes(x=EV1,y=EV2,col=lat),size=3, alpha=0.5) +
    scale_color_viridis(name="Latitude") + theme_bw() +
    ylab(paste("PC",j," (",round(pc.percent[j],3)," PVE)",sep="")) +
    xlab(paste("PC",i," (",round(pc.percent[i],3)," PVE)",sep=""))
  long <- ggplot(tab) + geom_point(aes(x=EV1,y=EV2,col=long),size=3, alpha=0.5) +
    scale_color_viridis(name="Longitude") + theme_bw() +
    ylab(paste("PC",j," (",round(pc.percent[j],3)," PVE)",sep="")) +
    xlab(paste("PC",i," (",round(pc.percent[i],3)," PVE)",sep=""))
  alt <- ggplot(tab) + geom_point(aes(x=EV1,y=EV2,col=alt),size=3,alpha=0.5) +
    scale_color_viridis(name="Altitude") + theme_bw() +
    ylab(paste("PC",j," (",round(pc.percent[j],3)," PVE)",sep="")) +
    xlab(paste("PC",i," (",round(pc.percent[i],3)," PVE)",sep=""))
  print(
    grid.arrange(lat, long, alt, nrow= 1)
  )
}
dev.off()

#Saving PCAs
tab <- data.frame(name = pca$sample.id,
                  PC1 = pca$eigenvect[,1],  
                  PC2 = pca$eigenvect[,2],    
                  PC3 = pca$eigenvect[,3],  
                  PC4 = pca$eigenvect[,4],
                  PC5 = pca$eigenvect[,5],  
                  PC6 = pca$eigenvect[,6],
                  PC7 = pca$eigenvect[,7],  
                  PC8 = pca$eigenvect[,8],
                  PC9 = pca$eigenvect[,9],  
                  PC10 = pca$eigenvect[,10],
                  PC11 = pca$eigenvect[,11],  
                  PC12 = pca$eigenvect[,12],
                  PC13 = pca$eigenvect[,13],  
                  PC14 = pca$eigenvect[,14],
                  PC15 = pca$eigenvect[,15],  
                  PC16 = pca$eigenvect[,16],
                  PC17 = pca$eigenvect[,17],  
                  PC18 = pca$eigenvect[,18],
                  PC19 = pca$eigenvect[,19],  
                  PC20 = pca$eigenvect[,20],
                  stringsAsFactors = FALSE)

write_tsv(tab, "Argophyllus.tranche90.snp.gwas.90.bi.ldr0p2.pca.txt")


#Plotting PCA loading vectors
eigen_corr <- snpgdsPCACorr(pca, genofile, eig.which=1:20, snp.id=snpset.id, num.thread=10)

snps.list <- eigen_corr$snp.id
eigen_corr_tibble <- tibble( snp.id = snps.list)

pdf("Argophyllus.tranche90.snp.gwas.90.bi.ldr0p2.pca.loadings.pdf")
for (i in 1:20){
  eigen_corr_tibble$correlation <- eigen_corr$snpcorr[i,]
  print(
    ggplot(eigen_corr_tibble, aes(x=snp.id, y=correlation)) + geom_hex() +
      theme_bw() + scale_fill_viridis() + ylab("PCA Loading") +
      xlab("SNP index") + ggtitle(paste("PC ",i,sep=""))
  )
}
dev.off()

###Getting the sites selected by LD filtering

chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
pos  <- read.gdsn(index.gdsn(genofile, "snp.position"))

all_snps <- tibble(chr = as.character(chr),pos=as.numeric(pos))
selected_snps <- all_snps[as.vector(snpset.id),]

selected_snps %>% mutate(chr2 = paste("HanXRQChr",chr,sep="")) %>%
  select(chr2, pos) %>% rename(chr = chr2) -> selected_snps

write_delim(selected_snps, "Argophyllus.tranche90.snp.gwas.90.bi.ldr0p2.txt")


### Getting IBD values

ibd.robust <- snpgdsIBDKING(genofile, snp.id=snpset.id, num.thread=10 )
ibd.dat <- snpgdsIBDSelection(ibd.robust)

pdf("Argophyllus.tranche90.snp.gwas.90.bi.ldr0p2.kinship.pdf")
ggplot(ibd.dat, aes(ID1, ID2 )) +
  geom_tile(aes(fill = kinship)) + scale_fill_viridis(name="KING kinship") +
  theme_bw()

ggplot(ibd.dat, aes(IBS0, kinship ),alpha=0.5) +
  geom_point() + ylab("KING kinship") + xlab("Proportion with 0 IBD") + 
  theme_bw()
dev.off()

write_delim(ibd.dat, "Argophyllus.tranche90.snp.gwas.90.bi.ldr0p2.kinship.txt")

##############################################
#Without inversions
##############################################

folder <- "/media/owens/Copper/wild_gwas_2018/argophyllus/"
gds.file <- "Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.v3noinv.ldfilter.gds"
vcf.file <- "Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.v3noinv.ldfilter.vcf.gz"
snpgdsVCF2GDS(paste(folder,vcf.file,sep="/"), paste(folder,gds.file,sep="/"), method="biallelic.only", ignore.chr.prefix = "Ha412HOChr")

genofile <- snpgdsOpen(paste(folder,gds.file,sep="/"))
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,method="r", num.thread=10)
snpset.id <- unlist(snpset)

pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=10, 
                 eigen.cnt = 0)

#Write eigenvectors and eigenvalues
write(pca$eigenval, "Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.v3noinv.ldfilter.ldr0p2.eigenvalues.txt",
      ncol=length(pca$eigenval))


write.table(pca$eigenvect, "Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.v3noinv.ldfilter.ldr0p2.eigenvectors.txt",
            col.names = F,row.names = F)



pc.percent <- pca$varprop*100
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T)
pop_loc <- read_tsv("pop_loc_allnum.txt")
pop_loc %>% rename(population = pop) %>% inner_join(.,labels) -> labels

snpgdsClose(genofile)


#Plotting PCAs
pdf("Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.v3noinv.ldfilter.ldr0p2.pca.pdf",height=6,width=14)
for (i in seq(1,20,2)){
  j <- i +1
  tab <- data.frame(name = pca$sample.id,
                    EV1 = pca$eigenvect[,i],  
                    EV2 = pca$eigenvect[,j],    
                    
                    stringsAsFactors = FALSE)
  inner_join(tab, labels) -> tab
  lat <- ggplot(tab) + geom_point(aes(x=EV1,y=EV2,col=lat),size=3, alpha=0.5) +
    scale_color_viridis(name="Latitude") + theme_bw() +
    ylab(paste("PC",j," (",round(pc.percent[j],3)," PVE)",sep="")) +
    xlab(paste("PC",i," (",round(pc.percent[i],3)," PVE)",sep=""))
  long <- ggplot(tab) + geom_point(aes(x=EV1,y=EV2,col=long),size=3, alpha=0.5) +
    scale_color_viridis(name="Longitude") + theme_bw() +
    ylab(paste("PC",j," (",round(pc.percent[j],3)," PVE)",sep="")) +
    xlab(paste("PC",i," (",round(pc.percent[i],3)," PVE)",sep=""))
  alt <- ggplot(tab) + geom_point(aes(x=EV1,y=EV2,col=alt),size=3,alpha=0.5) +
    scale_color_viridis(name="Altitude") + theme_bw() +
    ylab(paste("PC",j," (",round(pc.percent[j],3)," PVE)",sep="")) +
    xlab(paste("PC",i," (",round(pc.percent[i],3)," PVE)",sep=""))
  print(
    grid.arrange(lat, long, alt, nrow= 1)
  )
}
dev.off()

#Saving PCAs
tab <- data.frame(name = pca$sample.id,
                  PC1 = pca$eigenvect[,1],  
                  PC2 = pca$eigenvect[,2],    
                  PC3 = pca$eigenvect[,3],  
                  PC4 = pca$eigenvect[,4],
                  PC5 = pca$eigenvect[,5],  
                  PC6 = pca$eigenvect[,6],
                  PC7 = pca$eigenvect[,7],  
                  PC8 = pca$eigenvect[,8],
                  PC9 = pca$eigenvect[,9],  
                  PC10 = pca$eigenvect[,10],
                  PC11 = pca$eigenvect[,11],  
                  PC12 = pca$eigenvect[,12],
                  PC13 = pca$eigenvect[,13],  
                  PC14 = pca$eigenvect[,14],
                  PC15 = pca$eigenvect[,15],  
                  PC16 = pca$eigenvect[,16],
                  PC17 = pca$eigenvect[,17],  
                  PC18 = pca$eigenvect[,18],
                  PC19 = pca$eigenvect[,19],  
                  PC20 = pca$eigenvect[,20],
                  stringsAsFactors = FALSE)

write_tsv(tab, "Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.v3noinv.ldfilter.ldr0p2.pca.txt")


