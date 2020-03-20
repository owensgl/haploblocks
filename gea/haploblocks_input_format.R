#!/usr/bin/Rscript
## This Rscript conver 012 genotype file into Baypass allele count file

setwd("~/Documents/input_data/annuus_HA412/inversions")

data_1<-read.table(file = "Ha412HO_inv.dec11.H.annuus.inversions.genotypes.v1.IDs.attached.input.R.txt", header = FALSE)
colnames(data_1)<-c("ind_id","EV1","EV2","geno_count","genus","chrom","mds","pop_id","inv_id")

#data_ordered<-tail(data_1[order(data_1$pop_id),c(1:9)])

pop_code<-unique(data_1$pop_id)
inv_code<-unique(data_1$inv_id)
 

#### generate Baypass format input file ####
zero_let <- array (NA, c(length(inv_code),2*length(pop_code)))
L = length(pop_code)*2
locus_positions=(2*(unique(round((1:(L-2))/2)))+1)


for (i in 1:length(pop_code)){ 
    for (j in 1:length(inv_code)){
  
    #count <- count + 1  
    aa<-data_1[data_1$pop_id==pop_code[i],]
    bb<-aa[aa$inv_id==inv_code[j],]
    zero_let[j,locus_positions[i]]<- nrow(bb[bb$geno_count==0,])*2  + nrow(bb[bb$geno_count==1,])
    zero_let[j,locus_positions[i]+1]<-nrow(bb[bb$geno_count==2,])*2  + nrow(bb[bb$geno_count==1,])
    
    
    }
}   

write.table(zero_let,file = "Ha412HO_inv.dec11.H.annuus.inversions.genotypes.v1.IDs.attached.Allele.count.sorted.BAYPASS.input.txt", row.names= FALSE, col.names = FALSE)

####### generate population allele frequency ####

zero_let <- array (NA, c(2*length(inv_code),length(pop_code)))
L = length(inv_code)*2
locus_positions=(2*(unique(round((1:(L-2))/2)))+1)

    
for (i in 1:length(pop_code)){ 
    for (j in 1:length(inv_code)){
  
    #count <- count + 1  
    aa<-data_1[data_1$pop_id==pop_code[i],]
    bb<-aa[aa$inv_id==inv_code[j],]
    zero_let[locus_positions[j],i]<- (nrow(bb[bb$geno_count==0,])*2  + nrow(bb[bb$geno_count==1,]))/(nrow(bb)*2)
    zero_let[locus_positions[j]+1,i]<-(nrow(bb[bb$geno_count==2,])*2  + nrow(bb[bb$geno_count==1,]))/(nrow(bb)*2)
    
    
    }
}    

write.table(zero_let,file = "Ha412HO_inv.dec11.H.annuus.inversions.genotypes.v1.IDs.attached.sorted.BayEnv.allele.frequency.txt", row.names= FALSE, col.names = FALSE)
