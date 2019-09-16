library(tidyverse)

ft_coverage <- read_tsv("HaFT1.depth.txt",col_names = c("sample",	"gene",	"chr",	"start",
                                                        "end",	"length",	"raw_depth",	"raw_frac_mis"))

info <- read_tsv("sample_info_apr_2018.tsv") %>% rename(sample = name)

inv <- read_tsv("MDS_outliers/Ha412HO/argophyllus/Ha412HO_inv.jan09.pcasites.Ha412HOChr06.pos1.genotypes.txt")

pdf("Arg.Ha412HOChr06pos1.HaFT1.coverage.pdf",height=4,width=4)
ft_coverage %>%
  inner_join(.,info) %>%
  inner_join(.,inv) %>%
  filter(species == "Arg",!is.na(cluster_genotype)) %>%
  ggplot(.,aes(x=as.factor(triangle_genotype),y=1-raw_frac_mis)) + geom_boxplot() +
  theme_bw() + xlab("Ha412HOChr6.pos1") + ylab("Fraction of gene sequenced") 
dev.off()
