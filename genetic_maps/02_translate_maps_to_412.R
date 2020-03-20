###############################################################
### --- Translate all maps to Ha412 bp for use in syntR --- ###
###############################################################

library("dplyr")
library("stringr")
library("syntR")

### --- Pet pet and pet fal maps from Ostevik et al. 2019 --- ###

# define a function that will merge the two data types
convert2Ha412_1 <- function(oldpos, pos412) {
  left_join(oldpos, pos412, by = "map1_name") %>% 
    select(., map1_name, new_chr, new_pos, map2_name, map2_chr, map2_pos) %>%
    filter(., complete.cases(.)) %>% 
    mutate(new_chr = gsub("Ha412HOChr", "A", new_chr)) %>%
    filter(., new_chr %in% c("A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10", "A11", "A12", "A13", "A14", "A15", "A16", "A17")) %>%
    rename(., map1_chr = new_chr, map1_pos = new_pos) %>%
    return()
}

# neg
oldpos <- read.table("data/map_pet_fal_KO.txt", header = TRUE)
pos412 <- read.table("data/map_pet_fal_KO_Ha412_translate.txt", header = FALSE)
names(pos412) <- c("new_chr", "new_pos", "map1_name")
convert2Ha412_1(oldpos, pos412) %>% write.table(., file="data/412map_pet_fal.txt", quote = FALSE, row.names = FALSE)

# pet
oldpos <- read.table("data/map_pet_pet_KO.txt", header = TRUE)
pos412 <- read.table("data/map_pet_pet_KO_Ha412_translate.txt", header = FALSE)
names(pos412) <- c("new_chr", "new_pos", "map1_name")
convert2Ha412_1(oldpos, pos412) %>% write.table(., file="data/412map_pet_pet.txt", quote = FALSE, row.names = FALSE)

### --- ann maps from glycophate mapping project --- ###

# use same function as above
pos412 <- read.table("data/map_gly_Ha412_translate.txt", header = FALSE)
names(pos412) <- c("new_chr", "new_pos", "map1_name")
# gly1
oldpos <- read.table("data/map_ann_gly1_MT.txt", header = TRUE)
convert2Ha412_1(oldpos, pos412) %>% write.table(., file="data/412map_ann_gly1.txt", quote = FALSE, row.names = FALSE)
# gly2
oldpos <- read.table("data/map_ann_gly2_MT.txt", header = TRUE)
convert2Ha412_1(oldpos, pos412) %>% write.table(., file="data/412map_ann_gly2.txt", quote = FALSE, row.names = FALSE)
# gly3
oldpos <- read.table("data/map_ann_gly3_MT.txt", header = TRUE)
convert2Ha412_1(oldpos, pos412) %>% write.table(., file="data/412map_ann_gly3.txt", quote = FALSE, row.names = FALSE)


### --- Arg map from Barb et al. 2014 --- ###

convert2Ha412_2 <- function(oldpos, SFWpos) {
  left_join(SFWpos, oldpos, by = "marker") %>%
    mutate(., map1_name = marker, map1_chr = CHR, map1_pos = pos, map2_name = marker, map2_chr = chr, map2_pos = position) %>%
    select(map1_name, map1_chr, map1_pos, map2_name, map2_chr, map2_pos) %>%
    filter(., complete.cases(.)) %>%
    mutate(map1_chr = gsub("Ha412HOChr", "A", map1_chr)) %>%
    filter(., map1_chr %in% c("A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10", "A11", "A12", "A13", "A14", "A15", "A16", "A17")) %>%
    return()
}

# load SFW loci positions
SFWpos <- read.table("data/SFW_loci_412_filtered.txt", header=FALSE)
names(SFWpos) <- c("marker", "CHR", "pos")
SFWpos$marker <- gsub("SFW", "SFW0", SFWpos$marker)

# arg
oldpos <- read.table("data/map_arg_barb.txt", header = TRUE)
convert2Ha412_2(oldpos, SFWpos) %>% write.table(., file="data/412map_arg.txt", quote = FALSE, row.names = FALSE)


### --- Maps from Kai --- ###

convert2Ha412_3 <- function(map) {
  map %>% mutate(chr = recode(map$chrom, !!! setNames(paste0("L", 1:length(names(table(map$chrom)))), names(table(map$chrom))))) %>%
    mutate(CHR = str_split_fixed(map$markers, "_", 2)[,1], bp = str_split_fixed(map$markers, "_", 2)[,2]) %>% 
    mutate(., map1_name = markers, map1_chr = CHR, map1_pos = bp, map2_name = markers, map2_chr = chr, map2_pos = pos) %>%
    select(map1_name, map1_chr, map1_pos, map2_name, map2_chr, map2_pos) %>%
    filter(., complete.cases(.)) %>%
    mutate(map1_chr = gsub("Ha412HOChr", "A", map1_chr)) %>%
    filter(., map1_chr %in% c("A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10", "A11", "A12", "A13", "A14", "A15", "A16", "A17")) %>%
    return()
}
  
# dune
map <- convert2Ha412_3(read.table("data/map_pet_GSD_D_KH_new5.txt", header = TRUE))
map_list <- make_one_map(map)
plot_maps(map_df = map_list[[1]], map1_chrom_breaks = map_list[[2]], map2_chrom_breaks = map_list[[3]])
map_list <- make_one_map(map, flip_chrs = c("L1", "L2", "L6", "L7", "L11", "L14"))
plot_maps(map_df = map_list[[1]], map1_chrom_breaks = map_list[[2]], map2_chrom_breaks = map_list[[3]])
map_list[[1]] %>% mutate(chr = recode(.$map2_chr, !!! setNames(c("C05", "C01", "C17", "C13", "C14", "C08", "C11", "C10", "C03", "C04", 
                                                                 "C09", "C06.16", "C16.17", "C02A", "C02B", "C12.15", "C15", "C07"), 
                                                               paste0("L", 1:18)))) %>% 
  select(., map1_name, map1_chr, map1_pos, map2_name, chr, map2_pos) %>%
  rename(map2_chr = chr) %>% 
  write.table(., file="data/412map_pet_GSD_D.txt", quote = FALSE, row.names = FALSE)

# non-dune
map <- convert2Ha412_3(read.table("data/map_pet_GSD_N_KH.txt", header = TRUE))
map_list <- make_one_map(map)
plot_maps(map_df = map_list[[1]], map1_chrom_breaks = map_list[[2]], map2_chrom_breaks = map_list[[3]])
map_list <- make_one_map(map, flip_chrs = c("L1", "L2", "L3", "L4", "L8", "L10", "L11", "L16"))
plot_maps(map_df = map_list[[1]], map1_chrom_breaks = map_list[[2]], map2_chrom_breaks = map_list[[3]])
map_list[[1]] %>% mutate(chr = recode(.$map2_chr, !!! setNames(c("C01", "C07", "C11", "C15", "C08", "C02", "C06", "C09", "C04", "C17", 
                                                                 "C14", "C12.15", "C13", "C16.17", "C05", "C03", "C10"), 
                                                               paste0("L", 1:17)))) %>% 
  select(., map1_name, map1_chr, map1_pos, map2_name, chr, map2_pos) %>%
  rename(map2_chr = chr) %>% 
  write.table(., file="data/412map_pet_GSD_N.txt", quote = FALSE, row.names = FALSE)

# brook arg map2 (‘Southern map’ - parents = 1805.8 [early] and 1820.7 [late])
map <- read.table("data/map_arg_brook_map2.txt", header = TRUE)
map_list <- make_one_map(map)
plot_maps(map_df = map_list[[1]], map1_chrom_breaks = map_list[[2]], map2_chrom_breaks = map_list[[3]])
map <- read.table("data/map_arg_brook_map2.txt", header = TRUE) %>% filter(map2_chr < 16) %>% droplevels()
map_list <- make_one_map(map, flip_chrs = c("7", "9", "10", "11", "12", "13", "14"))
plot_maps(map_df = map_list[[1]], map1_chrom_breaks = map_list[[2]], map2_chrom_breaks = map_list[[3]])
map_list[[1]] %>% mutate(chr = recode(.$map2_chr, !!! setNames(c("C01", "C03", "C09", "C15.06", "C12.16", "C10", "C06.15", "C14",  
                                                                 "C17", "C11", "C02", "C07", "C04", "C05", "C13"), 
                                                               1:15))) %>% 
  select(., map1_name, map1_chr, map1_pos, map2_name, chr, map2_pos) %>%
  rename(map2_chr = chr) %>% 
  write.table(., file="data/412map_arg_brook2.txt", quote = FALSE, row.names = FALSE)


# brook arg map3 (Moyers et al. 2017, ‘Northern map’ -  parent = 1805.10 [early] and 1834.8 [late])
map <- read.table("data/map_arg_brook_map3.txt", header = TRUE)
map_list <- make_one_map(map)
plot_maps(map_df = map_list[[1]], map1_chrom_breaks = map_list[[2]], map2_chrom_breaks = map_list[[3]])
map_list <- make_one_map(map, flip_chrs = c("1", "2", "5", "6", "7", "8", "13", "14", "15", "16", "17"))
plot_maps(map_df = map_list[[1]], map1_chrom_breaks = map_list[[2]], map2_chrom_breaks = map_list[[3]])
map_list[[1]] %>% mutate(chr = recode(.$map2_chr, !!! setNames(c("C09", "C06.15", "C07.4", "C17", "C13", "C12.16", "C16.12", "C01",  
                                                                 "C04", "C11", "C03", "C10", "C06", "C02", "C08", "C14", "C05"), 
                                                               1:17))) %>% 
  select(., map1_name, map1_chr, map1_pos, map2_name, chr, map2_pos) %>%
  rename(map2_chr = chr) %>% 
  write.table(., file="data/412map_arg_brook3.txt", quote = FALSE, row.names = FALSE)
