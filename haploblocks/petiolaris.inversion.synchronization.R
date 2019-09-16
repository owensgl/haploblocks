#This script is for combining the pcasites for petiolaris inversions found in different sets of samples.
library(tidyverse)

map <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.synchronizePet.txt")

n <- 1

for (n in 1:nrow(map)){
  result <- NULL
  petpf <- NULL
  petpet <- NULL
  petfal <- NULL
  chr <- map$Chr[n]
  if (!is.na(map$PetPet[n])){
    petpet <- read_tsv(paste("MDS_outliers/Ha412HO/petpet/Ha412HO_inv.jan09.petpet.",chr,".",map$PetPet[n],".pcasites.txt.gz",sep=""))
  }
  if (!is.na(map$PetFal[n])){
    petfal <- read_tsv(paste("MDS_outliers/Ha412HO/petfal/Ha412HO_inv.jan09.petfal.",chr,".",map$PetFal[n],".pcasites.txt.gz",sep=""))
  }
  if (!is.na(map$PetPF[n])){
    petpf <- read_tsv(paste("MDS_outliers/Ha412HO/petiolaris/Ha412HO_inv.jan09.petiolaris.",chr,".",map$PetPF[n],".pcasites.txt.gz",sep=""))
  }

  if (!is.null(petpf)){
    result <- petpf
  }
  if (!is.null(result) & !is.null(petpet)){
    result <- result %>% inner_join(petpet %>% select(xrqchr,xrqpos))
  }else if (!is.null(petpet)){
    result <- petpet
  }
  
  if (!is.null(result) & !is.null(petfal)){
    result <- result %>% inner_join(petfal %>% select(xrqchr,xrqpos))
  }else if (!is.null(petfal)){
    result <- petfal
  }
  write_tsv(result,paste("MDS_outliers/Ha412HO/petiolaris_syn/Ha412HO_inv.jan09.petiolaris_syn.",chr,".",map$synchronized[n],".pcasites.txt.gz",sep=""))
  print(c(map$Chr[n],nrow(petpf),nrow(petpet),nrow(petfal), nrow(result)))
  
}


