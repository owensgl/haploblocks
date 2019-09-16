#This is for making a figure of inversion locations
library(tidyverse)
library(statebins)
chr_sizes <- read_tsv("Ha412HO.chrlengths.txt")

inv_locations <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt")

inv_locations %>%
  separate(chr, c("ref","chr_n"), "Chr") %>%
  mutate(chr_n = as.numeric(chr_n)) %>%
  filter(species != "Niveus") %>%
  mutate(species_vert = case_when(species == "Annuus" ~ 0.66,
                                  species == "Argophyllus" ~ 0.33,
                                  species == "Petiolaris" ~ 0),
         full_name = case_when(species == "Annuus" ~ "H. annuus",
                               species == "Argophyllus" ~ "H. argophyllus",
                               species == "Petiolaris" ~ "H. petiolaris")) -> inv_locations_formatted
  

  
pdf("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.pdf",height=5,width=5)
chr_sizes %>%
  separate(chr, c("ref","chr_n"), "Chr") %>%
  mutate(chr_n = as.numeric(chr_n)) %>%
  ggplot(.) + 
  geom_rect(data=inv_locations_formatted,aes(xmin=start/1000000,xmax=end/1000000,ymin=chr_n+species_vert-0.5,ymax=chr_n+species_vert+0.33-0.5,fill=full_name),color=NA) +
  geom_rect(aes(xmin=start/1000000,xmax=end/1000000,ymin=chr_n-0.5,ymax=chr_n+0.5),color="black",fill=NA) +
  scale_fill_manual(values =c("#FFB31C","#3DA047","#447499"),name="Species",
                    labels = c(expression(italic("H. annuus")),expression(italic("H. argophyllus")),expression(italic("H. petiolaris")))) +
  scale_y_reverse(breaks=c(1:17)) +
  xlab("MBp") + ylab("Ha412HO Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        legend.position="bottom",
        axis.text.y = element_text(margin=margin(0,-15,0,0)),
        axis.ticks.y = element_blank())
dev.off()
  



