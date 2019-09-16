library(tidyverse)
library(lostruct)
library(Matrix)
library(colorspace)
library(RColorBrewer)
library(ggmap)
library(scatterpie)
library(SNPRelate)
library(gridExtra)
library(ggExtra)
library(grid)
library(zoo)
for (chr in sprintf("%02d", seq(1,17))){
  genome_prefix <- "Ha412HOChr"
  data_directory <- "/media/owens/Copper/wild_gwas_2018/argophyllus/"
  data_name <- paste("Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO",sep="")
  
  labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T)
  pop_loc <- read_tsv("pop_loc_allnum.txt")
  pop_loc %>% rename(population = pop) %>% inner_join(.,labels) -> labels
  
  #Prepare map data
  usa <- map_data('state')
  states <- map_data("state")
  target_state <- map_data('state')
  lat_range <- c(25, 50)
  long_range <- c(-125,-93)
  pie_size <- 0.4
  
  window_size <- 100
  k_kept <- 40
  max_distance_between_outliers <- 100
  
  bcf.file <- paste(data_directory, "/", data_name,".Chr",chr, ".bcf",sep="")
  samples <- read_tsv(paste(data_directory, "/", data_name,".samplelist.txt",sep=""),col_names = F)
  colnames(samples) <- c("sample")
  sites <- vcf_positions(bcf.file)
  win.fn.snp <- vcf_windower(bcf.file, size=window_size, type="snp", sites=sites) 
  system.time( snp.pca <- eigen_windows(win.fn.snp,k=2, mc.cores=10) )
  system.time( pcdist <- pc_dist( snp.pca ) )
  
  pcdist_na <- which(is.na(pcdist), TRUE)
  
  
  na.inds <- is.na(pcdist[,1]) 
  if (sum(na.inds) == length(na.inds)){
    na.inds <- is.na(pcdist[,2]) 
  }
  mds <- cmdscale( pcdist[!na.inds,!na.inds], eig=TRUE, k=k_kept )
  
  mds.coords <- mds$points
  colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))
  win.regions <- region(win.fn.snp)()
  win.regions$n <- 1:nrow(win.regions)
  win.regions <- win.regions[!na.inds,]
  win.regions %>% mutate(mid = (start + end) / 2) ->  win.regions
  #Add the columns for all the MDS coordinates
  for (k in 1:k_kept){
    str_pad(k, 2, pad = "0")
    
    name = paste("mds",str_pad(k, 2, pad = "0"),sep="")
    win.regions$tmp <- "NA"
    win.regions <- win.regions %>% rename(!!name := tmp) 
  }
  
  #Add the MDS coordinates to each window.
  for (i in 1:k_kept){
    j = i + 5
    win.regions[,j] <- mds.coords[,i]
  }
  
  #Make plots of all the MDS to visualize patterns
  pdf(paste(data_name, ".Chr", chr, ".",window_size,".MDSplots.pdf",sep=""), height=16,width=25)
  print(
    win.regions %>%
      gather(., mds, value, colnames(win.regions)[6:(ncol(win.regions)-30)]) %>% 
      ggplot(.,aes(x=mid,y=value)) + geom_point() + facet_grid(mds~.,scales = "free") +
      theme_bw()
  )
  
  dev.off()
  saveRDS(win.regions, file = paste(data_name, ".Chr", chr, ".",window_size,".windows.rds",sep=""))
}