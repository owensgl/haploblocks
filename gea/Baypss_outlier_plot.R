#!/usr/bin/env Rscript

library(qqman)

setwd("/data/home/shaghayegh/H_annuus_2018_HA412/baypass/maf_0.03/downsteram_analyses")
all_good_chrom<-read.table(file = "var_out_annuus_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.txt",header = FALSE)
colnames(all_good_chrom)<-c("snp_id","CHR","pos","latitude_e","longitude_e","elevation_e","MAT_e","MWMT_e","MCMT_e","TD_e","MAP_e","MSP_e","AHM_e","SHM_e","DD_0_e","DD5_e","DD_18_e","DD18_e","NFFD_e","bFFP_e","eFFP_e","FFP_e","PAS_e","EMT_e","EXT_e","Eref_e","CMD_e","MAR_e","RH_e")


test_type_env<-array ("0",ncol (all_good_chrom))
env_type<-grep ("_e$",colnames(all_good_chrom))
test_type_env[min(env_type):(min(env_type)+25)]<-"envir"
all_good_chrom$id<- paste (all_good_chrom[,2],all_good_chrom[,3], sep = "__")


#### Top candidate windows and BF threshold ####
out_res <- NULL
test_names_env <- colnames (all_good_chrom)[test_type_env == "envir"]
the_i <- which (test_type_env == "envir") 
bf <- 10


## extract outliers (here SNPs with BFs more than significance threshold based on Baypass paper )
for (i in 1:length (the_i)){

	outliers_snp <- (which(all_good_chrom[,the_i[i]] >= bf))
	snps <- (which(all_good_chrom[,the_i[i]] < 10000))

    outlier_count<-length(outliers_snp)
    snps_count <- length(snps)
    
    ## get the percentage
    outlier_percent<-((outlier_count/(snps_count))*100)

   sub_good <- data.frame (snps_count,outlier_count,outlier_percent)
   sub_good$test_name <- test_names_env[i]		
   out_res <- rbind (out_res,sub_good)

}

write.table(out_res, file = "outlire_percentage_annuus_2018_HA412_baypass_STD_IS2_bf10_whole_genome.txt", col.names = TRUE, row.names = FALSE)


for(i in 1:length(test_names_env)) {
pdf(file = paste("manhattan_plots_baypass_STD_IS2_annuus_2018_HA412_whole_genome",test_names_env[i], ".pdf", sep=""), width = 24)
manhattan(all_good_chrom,bp ="pos",chr = "CHR", p = test_names_env[i], snp = "snp_id", logp = FALSE, ylim = c(0,60), genomewideline = 10, suggestiveline = FALSE ,xlab = "chromosome",
ylab = "BFis (in dB)", cex.axis=1.5, cex.lab=1.4,col = c("blue4","orange3"),main = test_names_env[i])     
dev.off()
}


chromosome<-unique(all_good_chrom$CHR)
for(i in 1:length(test_names_env)){

for (j in 1:length(chromosome)){
pdf(file = paste("manhattan_plots_baypass_STD_IS2_annuus_HA412_chrom",chromosome[j],"_2018_all_",test_names_env[i], ".pdf", sep=""), width = 24)
manhattan(subset(all_good_chrom,CHR == chromosome[j]), bp ="pos", p = test_names_env[i], snp = "snp_id", logp = FALSE, ylim = c(0,60), genomewideline = 10, suggestiveline = FALSE ,xlab = "chromosome",
ylab = "BFis (in dB)", cex.axis=1.5, cex.lab=1.4,col = c("blue4","orange3"),main = paste(test_names_env[i],"_",chromosome[j]))
dev.off()
}
}
