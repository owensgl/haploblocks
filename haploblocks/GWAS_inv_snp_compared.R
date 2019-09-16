#This compares GWAS inversion association to SNPs
library(tidyverse)

pvalue_cutoff <- 0.001
min_freq <- 0.03
directory <- "/media/owens/Copper/wild_gwas/gwas"
species_abbreviation_list <- c("annuus","argophyllus","petpet","petfal")
species_capital_list <-c("Annuus", "Argophyllus","Petiolaris","Petiolaris")
tag_list <- c("gwas", "gwas","petpet","petfal")



gwas_scores <- tibble(species=character(),threshold=character(),
                      inv_positive=numeric(),inv_count=numeric(),
                      snp_positive=numeric(),snp_count=numeric(),
                      pvalue=numeric())
for (n in 1:4){
  
  
  
  
  species_abbreviation <- species_abbreviation_list[n]
  species_capital <- species_capital_list[n]
  tag <- tag_list[n]
  
  ##Get list of all traits
  traits <- tibble(trait=list.files(paste0("gwas/",species_abbreviation),include.dirs = F) ) %>%
    filter(!grepl("znorm",trait)) %>%
    separate(trait,c("trait","spacer"),"\\.") %>%
    select(-spacer)
  
  #Get the list of allele frequencies
  freqs <- read_tsv(paste0(directory,"/",species_abbreviation,"/",
                           species_capital,".tranche90.snp.",tag,".90.bi.remappedHa412HO.beagle.freq"),
                    col_names=c("chr","pos","freq")) %>%
    filter(freq > min_freq, freq < (1-min_freq)) %>%
    select(-freq)
  
  
  
  #Get the list off full sites and remove sites in LNH
  all_snps <- read_tsv(paste0(directory,"/",species_abbreviation,"/results/",
                              species_capital,".tranche90.snp.",tag,".90.bi.remappedHa412HO.",
                              traits$trait[1],".beagle.v3noinv.ldfilter.ps.gz"),
                       col_names=c("id","beta","sd","pvalue")) %>%
    select(id) %>% separate(id,c("chr","pos"),"_") %>%
    mutate(pos=as.numeric(pos)) %>%
    semi_join(.,freqs) %>%
    mutate(count=0)
  
  inversion_regions <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt")
  chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")
  
  inversions <- chr_lengths %>% 
    select(chr,end) %>%
    rename(chr_len = end) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(inversion_regions, ., by=c("chr"="chr")) %>%
    # Add a cumulative position of each SNP
    arrange(chr, start) %>%
    mutate( startcum=start+tot,endcum=end+tot) %>%
    filter(species == species_capital) 
  
  all_snps  %>%
    mutate(pos=as.numeric(pos))-> all_sites
  
  all_sites %>%
    mutate(window = floor(pos/1000000),
           window_start = window*1000000,
           window_end = (window+1)*1000000) %>%
    select(chr,window,window_start,window_end) %>%
    unique() -> all_windows
  
  windows_in_inversions <- tibble(chr=character(), window=character(), 
                                  window_start=numeric(),window_end=numeric())
  for (i in 1:nrow(inversions)){
    chosen_chr = gsub("Ha412HOChr","",inversions$chr[i],)
    chosen_start = inversions$start[i]
    chosen_end = inversions$end[i]
    all_windows %>%
      filter(chr == chosen_chr,window_end > chosen_start, window_start < chosen_start) -> left_side
    all_windows %>%
      filter(chr == chosen_chr,window_end > chosen_end, window_start < chosen_end) -> right_side 
    all_windows %>%
      filter(chr == chosen_chr,window_end < chosen_end, window_start > chosen_start) -> middle
    all_windows %>%
      filter(chr == chosen_chr,window_end > chosen_end, window_start < chosen_start) -> surrounding
    
    
    windows_in_inversions <- rbind(windows_in_inversions, left_side, right_side, middle, surrounding)
    
  }
  
  windows_outside_inversions <- all_windows %>%
    anti_join(., windows_in_inversions)
  
  
  all_snps %>%
    mutate(pos=as.numeric(pos)) %>%
    mutate(window = floor(pos/1000000),
           window_start = window*1000000,
           window_end = (window+1)*1000000) %>%
    semi_join(.,windows_outside_inversions) %>%
    select(-window,-window_start,-window_end)-> all_snps
  
  for (x in 1:nrow(traits)){
    gwas_snps <- read_tsv(paste0(directory,"/",species_abbreviation,"/results/",
                                 species_capital,".tranche90.snp.",tag,".90.bi.remappedHa412HO.",
                                 traits$trait[x],".beagle.v3noinv.ldfilter.ps.gz"),
                          col_names=c("id","beta","sd","pvalue"))
    
    
    
    
    gwas_snps %>%
      separate(id,c("chr","pos"),"_") %>%
      mutate(pos=as.numeric(pos)) %>%
      mutate(window = floor(pos/1000000),
             window_start = window*1000000,
             window_end = (window+1)*1000000) %>%
      semi_join(.,windows_outside_inversions)  %>%
      select(-window,-window_start,-window_end) %>%
      filter(pvalue < pvalue_cutoff) %>%
      select(chr,pos) %>%
      semi_join(.,freqs) %>%
      mutate(newcount=1)-> hits
    
    all_snps <- full_join(all_snps,hits) %>% mutate(newcount = replace_na(newcount, 0)) %>%
      mutate(summed_count = count+newcount) %>%
      select(-count,-newcount) %>% rename(count=summed_count)
  }
  
  
  
  
  #Pull out associations for inversions
  
  all_invs <- read_tsv(paste0(directory,"/",species_abbreviation,"/results/",
                              "Ha412HO_inv.v3.pcasites.",species_abbreviation,
                              ".shuffled.",
                              traits$trait[1],".ps.gz"),
                       col_names=c("id","beta","sd","pvalue")) %>%
    select(id) %>% separate(id,c("chr","pos"),"_") %>%
    mutate(count=0)
  
  for (x in 1:nrow(traits)){
    gwas_invs <- read_tsv(paste0(directory,"/",species_abbreviation,"/results/",
                                 "Ha412HO_inv.v3.pcasites.",species_abbreviation,
                                 ".shuffled.",
                                 traits$trait[x],".ps.gz"),
                          col_names=c("id","beta","sd","pvalue"))
    gwas_invs %>%
      filter(pvalue < pvalue_cutoff) %>%
      select(id) %>%
      separate(id,c("chr","pos"),"_") %>%
      mutate(newcount=1)-> hits
    
    all_invs <- full_join(all_invs,hits) %>% mutate(newcount = replace_na(newcount, 0)) %>%
      mutate(summed_count = count+newcount) %>%
      select(-count,-newcount) %>% rename(count=summed_count)
  }
  ##Remove inversions below a minimum allele frequency
  if (species_capital == "Petiolaris"){
    inv_frequencies <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.petinvfreq.txt") %>%
      filter(tolower(species) == species_abbreviation) %>%
      mutate(chr = gsub("Ha412HOChr","",chr),
             pos=gsub("syn","",mds)) %>%
      select(-mds,-species)
    inner_join(all_invs, inv_frequencies) %>%
      filter(freq > min_freq,freq < (1-min_freq)) %>%
      select(-freq) -> all_invs
  }
  
  
  
  ##Test if one or more phenotypes are associated
  inv_test <- all_invs %>%
    group_by(count) %>%
    summarize(total=n()) %>%
    mutate(any_sig = case_when( count >= 1 ~ "Yes",
                                TRUE ~ "No")) %>%
    group_by(any_sig) %>%
    summarize(count=sum(total))
  
  snp_test <- all_snps %>%
    group_by(count) %>%
    summarize(total=n()) %>%
    mutate(any_sig = case_when( count >= 1 ~ "Yes",
                                TRUE ~ "No")) %>%
    group_by(any_sig) %>%
    summarize(count=sum(total))
  
  test_result <- prop.test(c(inv_test$count[which(inv_test$any_sig == "Yes")], snp_test$count[which(snp_test$any_sig == "Yes")]),
                           c(sum(inv_test$count),sum(snp_test$count)))
  
  tmp.tibble <- tibble(species=species_abbreviation,threshold="1+",
                       inv_positive=inv_test$count[which(inv_test$any_sig == "Yes")],
                       inv_count=sum(inv_test$count),
                       snp_positive=snp_test$count[which(snp_test$any_sig == "Yes")],
                       snp_count=sum(snp_test$count),pvalue=test_result$p.value)
  gwas_scores <- rbind(gwas_scores,tmp.tibble)
  
  ####Test if 2 or more phenotypes are associated
  inv_test <- all_invs %>%
    group_by(count) %>%
    summarize(total=n()) %>%
    mutate(any_sig = case_when( count >= 2 ~ "Yes",
                                TRUE ~ "No")) %>%
    group_by(any_sig) %>%
    summarize(count=sum(total))
  
  snp_test <- all_snps %>%
    group_by(count) %>%
    summarize(total=n()) %>%
    mutate(any_sig = case_when( count >= 2 ~ "Yes",
                                TRUE ~ "No")) %>%
    group_by(any_sig) %>%
    summarize(count=sum(total))
  
  test_result <- prop.test(c(inv_test$count[which(inv_test$any_sig == "Yes")], snp_test$count[which(snp_test$any_sig == "Yes")]),
                           c(sum(inv_test$count),sum(snp_test$count)))
  
  tmp.tibble <- tibble(species=species_abbreviation,threshold="2+",
                       inv_positive=inv_test$count[which(inv_test$any_sig == "Yes")],
                       inv_count=sum(inv_test$count),
                       snp_positive=snp_test$count[which(snp_test$any_sig == "Yes")],
                       snp_count=sum(snp_test$count),pvalue=test_result$p.value)
  gwas_scores <- rbind(gwas_scores,tmp.tibble)
}

write_tsv(gwas_scores,"gwas/Ha412HO_inv.v3.gwas.LNHenrichment.txt")
gwas_scores %>%
  mutate(percent_inv = inv_positive/inv_count,
         percent_snp = snp_positive/snp_count) %>%
  gather(type,value,percent_inv:percent_snp) %>%
  mutate(species = case_when(species=="annuus" ~ "H. annuus",
                             species=="argophyllus"~"H. argophyllus",
                             species=="petpet"~"H. petiolaris petiolaris",
                             species=="petfal"~"H. petiolaris fallax")) %>%
  mutate(threshold=case_when(threshold=="1+"~"1+ variables",
                             threshold=="2+"~"2+ variables")) %>%
  mutate(sig_dot=case_when(pvalue < 1e-20 ~ "**",
                           pvalue < 0.05 ~ "*",
                           TRUE ~ " ")) %>%
  group_by(species,threshold) %>%
  mutate(max_prop = max(value))%>%
  ggplot(.,aes(x=species,y=value,fill=type)) + geom_bar(stat="identity",position="dodge") +
  geom_text(aes(x=species,y=max_prop+0.02,label=sig_dot)) +
  facet_wrap(~threshold) +
  ylab("Proportion associated") + xlab("Species") +
  theme_bw() + 
  scale_fill_manual(values=c("grey","black"),name="Variant Type",
                    labels=c("LNH","SNP")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("gwas/Ha412HO_inv.v3.gwas.LNHenrichment.pdf")
