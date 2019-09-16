# README:
# This uses the following genetic map files:
# CR29xHA468.csv
# HA412-HOxANN1238.csv
# HA412-HOxRHA415.csv
# HA89xRHA464.csv
# NMS373xHopi.csv
# Nov22k22.consolidated_split.narrow.hardmask_split.best_cM_Ha412HOv2.0-20181130_sorted.csv
# NuSdBxHA464.csv
# RHA280xRHA801.csv
# 
# It creates Ha412HOv2.0-20181130.Nov22k22.geneticmap.csv and Ha412HOv2.0-20181130.Nov22k22.geneticmap.txt
# The script divide_up_cm_map.pl is used to create Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt


library(tidyverse)
library(RColorBrewer)
directory <- "/media/owens/Copper/wild_gwas_2018/genetic_maps"
chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")
chr_lengths$LG <-  as.numeric(gsub("Ha412HOChr","",chr_lengths$chr))
files <- list.files(directory)

maps <- tibble(chr=character(),pos=numeric(),LG=character(),cM=numeric(),
               map=character())
for (i in 1:length(files)){
  data <- read_csv(paste0(directory,"/",files[i]),col_names = c("chr","pos","LG","cM"),skip=1)
  data$map <- gsub(".csv","",files[i])
  data$cM <- as.numeric(data$cM)
  data$LG <- as.numeric(data$LG)
  maps <- rbind(maps,data)
}
#Flip some chr 12 for 3 samples
maps %>%
  mutate(chr_lg = as.numeric(gsub("Ha412HOChr","",chr))) %>%
  filter(chr_lg == LG) %>%
  group_by(map,chr) %>%
  mutate(max_cM = max(cM)) %>%
  ungroup() %>%
  mutate(cM = case_when(chr == "Ha412HOChr12" & map == "CR29xHA468"~ abs(cM - max_cM),
                        chr == "Ha412HOChr12" & map == "HA89xRHA464"~ abs(cM - max_cM),
                        chr == "Ha412HOChr12" & map == "NuSdBxHA464"~ abs(cM - max_cM),
                        chr == "Ha412HOChr15" & map == "CR29xHA468"~ abs(cM - max_cM),
                        chr == "Ha412HOChr15" & map == "HA89xRHA464"~ abs(cM - max_cM),
                        chr == "Ha412HOChr15" & map == "NuSdBxHA464"~ abs(cM - max_cM),
                        TRUE ~ cM)) -> maps

  

inversions <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.annuus.bed",
                       col_names = c("chr","start","stop")) %>%
  mutate(LG = as.numeric(gsub("Ha412HOChr","",chr)))


pdf("HA412HO_inv.v3.geneticmaps.pdf",height=12,width=16)
maps %>%
  filter(map !="Ha412HOv2.0-20181130.Nov22k22.geneticmap") %>%
  filter(map !="Nov22k22.consolidated_split.narrow.hardmask_split.best_cM_Ha412HOv2.0-20181130_sorted") %>%
  ggplot(.,) + geom_line(aes(x=pos/1000000,y=cM,color=map),size=1,
                         alpha=0.5) +
  scale_color_brewer(palette = "Set1") +
  geom_segment(data=inversions,aes(x=start/1000000,xend=stop/1000000,y=-5,yend=-5),
               size=2) +
  geom_line(data=maps %>% 
              filter(map =="Ha412HOv2.0-20181130.Nov22k22.geneticmap"),
            aes(x=pos/1000000,y=cM),color="Black",size=1.4) +
  theme_bw() + facet_wrap(~LG) +
  xlab("MB")

maps %>%
  filter(map=="HA412-HOxANN1238") %>%
  ggplot(.,) + geom_line(aes(x=pos/1000000,y=cM),size=1,color=brewer.pal(3,"Set1")[2]) +
  geom_segment(data=inversions,aes(x=start/1000000,xend=stop/1000000,y=-5,yend=-5),size=2) +
  theme_bw() + facet_wrap(~LG) +
  xlab("MB") + ggtitle("HA412-HOxANN1238")

dev.off()

###Using smooth spline to remove outliers

interpolated_map <- tibble(chr=character(),pos=numeric(),cM=numeric() )

pdf("Ha412HOv2.0-20181130.Nov22k22.geneticmap.pdf")
for (chr_n in seq(17)){
  chr_chosen <- chr_lengths$chr[which(chr_lengths$LG == chr_n)]
  
  positions <- maps %>% filter(map=="Nov22k22.consolidated_split.narrow.hardmask_split.best_cM_Ha412HOv2.0-20181130_sorted") %>%
    filter(LG == chr_n) %>%
    arrange(pos) %>%
    pull(pos) 
  
  cM <- maps %>% filter(map=="Nov22k22.consolidated_split.narrow.hardmask_split.best_cM_Ha412HOv2.0-20181130_sorted") %>%
    filter(LG == chr_n) %>%
    arrange(pos) %>%
    pull(cM) 
  #Flip the start of Chr12 because of likely assembly issue
  if (chr_n == 12){
    positions <- maps %>% filter(map=="Nov22k22.consolidated_split.narrow.hardmask_split.best_cM_Ha412HOv2.0-20181130_sorted") %>%
      filter(LG == chr_n) %>%
      arrange(pos) %>%
      mutate(pos = case_when(pos <= 10213302 ~ abs(pos - 10213301),
                                  TRUE ~ pos)) %>%
      pull(pos)
  }
  original_tibble = tibble(positions=positions,cM=cM)
  spline <- smooth.spline(x=positions,y=cM,spar=0.5)
  
  filtered_positions <- positions[which(abs(residuals(spline)) < 1)]
  filtered_cM <- cM[which(abs(residuals(spline)) < 1)]
  spline_filtered <- smooth.spline(x=filtered_positions,y=filtered_cM,spar=0.5)
  spline_tibble = tibble(pos=spline_filtered$x,cM=spline_filtered$y)
  
  
  for (i in c(seq(1, chr_lengths$end[which(chr_lengths$LG == chr_n)], 1000000),chr_lengths$end[which(chr_lengths$LG == chr_n)])){
    cM <- round(predict(spline_filtered, i)$y,2)
    if (cM < 0){
      cM = 0.01
    }
    if (i > 1){
      if (tail(interpolated_map,1)$cM > cM){
        cM <- tail(interpolated_map,1)$cM
      }
    }
    tmp <- tibble(chr=chr_chosen,pos=i, cM=cM)
    interpolated_map <- rbind(interpolated_map,tmp)
  }

  print(
    ggplot(spline_tibble,aes(x=pos/1000000,y=cM)) + geom_line() +
      theme_bw() + geom_line(data=interpolated_map %>% filter(chr == chr_chosen),
                             aes(x=pos/1000000,y=cM),color="red",alpha=0.8) +
      ggtitle(chr_chosen) + xlab("MB") +
      geom_point(data=original_tibble,aes(x=positions/1000000,y=cM),alpha=0.05)
  )
}
dev.off()
interpolated_map$LG <- gsub("Ha412HOChr","",interpolated_map$chr)
interpolated_map %>% select(chr,pos,LG, cM) -> interpolated_map
write_tsv(interpolated_map,"Ha412HOv2.0-20181130.Nov22k22.geneticmap.txt")

write_csv(interpolated_map,"/media/owens/Copper/wild_gwas_2018/genetic_maps/Ha412HOv2.0-20181130.Nov22k22.geneticmap.csv")


