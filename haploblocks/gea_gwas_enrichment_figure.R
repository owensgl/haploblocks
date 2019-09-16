library(tidyverse)


gwas_scores <- read_tsv("gwas/Ha412HO_inv.v3.gwas.LNHenrichment.txt") %>%
  mutate(data_type = "Trait GWA")
gea_scores <- read_tsv("gea/var_out_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.LNHenrichment.txt") %>%
  mutate(data_type = "Environmental GEA")

pdf("gwas/Ha412HO_inv.v3.gwas_gea.enrichment.pdf")
rbind(gwas_scores, gea_scores) %>%
  filter(threshold == "1" | threshold == "1+") %>%
  select(-threshold) %>%
  mutate(percent_inv = inv_positive/inv_count,
         percent_snp = snp_positive/snp_count) %>%
  gather(type,value,percent_inv:percent_snp) %>%
  mutate(species = case_when(species=="annuus" ~ "H. annuus",
                             species=="argophyllus"~"H. argophyllus",
                             species=="petpet"~"H. petiolaris petiolaris",
                             species=="petfal"~"H. petiolaris fallax")) %>%
  mutate(sig_dot=case_when(pvalue < 0.0005 ~ "**",
                           pvalue < 0.05 ~ "*",
                           TRUE ~ " ")) %>%
  group_by(species,data_type) %>%
  mutate(max_prop = max(value))%>%
  ggplot(.,aes(x=species,y=value,fill=type)) + geom_bar(stat="identity",position="dodge") +
  geom_text(aes(x=species,y=max_prop+0.02,label=sig_dot)) +
  facet_wrap(~data_type) +
  ylab("Proportion associated") + xlab("Species") +
  theme_minimal() + 
  scale_fill_manual(values=c("grey","black"),name="Variant Type",
                    labels=c("SV","SNP")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()
