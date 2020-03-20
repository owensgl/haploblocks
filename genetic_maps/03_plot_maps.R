########################################################
### --- Plot LNH regions on sets of genetic maps --- ###
########################################################

library("dplyr")
library("scales")
library("syntR")

# make a function to add SV locations to data
add_SVs <- function(map, intervals){
  positions <- vector()
  for (i in 1:nrow(intervals)) {
    positions <- c(positions, which(map$map1_pos>intervals$start[i] & map$map1_pos<intervals$end[i] & as.character(map$map1_chr)==as.character(intervals$chr[i]))
    )
  }
  map$color <- 1
  map$color[positions] <- 0
  map
}

# load intervals 
pet_intervals <- read.table("data/SV_locations_pet.txt", header = TRUE)
arg_intervals <- read.table("data/SV_locations_arg.txt", header = TRUE)
ann_intervals <- read.table("data/SV_locations_ann.txt", header = TRUE)

# load maps
pet1_map <- read.table("data/412map_pet_pet.txt", header = TRUE) %>% add_SVs(., pet_intervals)
pet2_map <- read.table("data/412map_pet_fal.txt", header = TRUE) %>% add_SVs(., pet_intervals)
pet3_map <- read.table("data/412map_pet_GSD_D.txt", header = TRUE) %>% add_SVs(., pet_intervals)
pet4_map <- read.table("data/412map_pet_GSD_N.txt", header = TRUE) %>% add_SVs(., pet_intervals)
arg1_map <- read.table("data/412map_arg.txt", header = TRUE) %>% add_SVs(., arg_intervals)
arg2_map <- read.table("data/412map_arg_brook3.txt", header = TRUE) %>% add_SVs(., arg_intervals)
arg3_map <- read.table("data/412map_arg_brook2.txt", header = TRUE) %>% add_SVs(., arg_intervals)
ann1_map <- read.table("data/412map_ann_gly1.txt", header = TRUE) %>% add_SVs(., ann_intervals)
ann2_map <- read.table("data/412map_ann_gly2.txt", header = TRUE) %>% add_SVs(., ann_intervals)
ann3_map <- read.table("data/412map_ann_gly3.txt", header = TRUE) %>% add_SVs(., ann_intervals)


map_list <- make_one_map(ann1_map)
plot_maps(map_df = map_list[[1]], map1_chrom_breaks = map_list[[2]], map2_chrom_breaks = map_list[[3]])


# plot a specific map
temp <- pet1_map %>% filter(map1_chr == "A17", map2_chr == "P16.17" | map2_chr == "P17")
temp <- pet1_map %>% filter(map1_chr == "A05", map2_chr == "P05")
plot(temp$map2_pos ~ temp$map1_pos, col = (temp$color + 1), 
     xlab = NA, ylab = NA, pch = 16, xlim = NULL, ylim = NULL, las = 2)

######################################
### --- Final petiolaris plots --- ###
######################################

layout(matrix(c(1,1,1,1,1,2,3,3,3,3,3,4,5,5,5,5,5,6,7,7,7,7,7,8), 4, 6, byrow = TRUE))
par(mai=c(0.2,0.2,0.2,0.2), mar=c(2, 4, 1.5, 1.5))
pet_map_list <- list(pet1_map, pet2_map, pet3_map, pet4_map)
# inversion, fallax, inversion, petiolaris #"#645084"
color_palette <- c("#953071", "#447499", "#953071", "#4BB7B7")
col_list <- c(3,1,1,1)

# pet01.01
chr_list <- c("P01", "N01", "C01", "C01")
y_list1 <- c(0,0,0,0)
y_list2 <- c(80,70,30,100)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A01", map2_chr == chr_list[i], map1_pos<110000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,110000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet01.01")
  polygon(c(180352,180352,12609141,12609141), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(90815284,90815284,98571815,98571815), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
  }

# pet05.01 
chr_list <- c("P05", "N05", "C05", "C05")
y_list1 <- c(50,80,100,50)
y_list2 <- c(100,130,170,100)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A05", map2_chr == chr_list[i], map1_pos>140000000, map1_pos<187000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(140000000,187000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet05.01")
  polygon(c(180352,180352,12609141,12609141), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(154541712,154541712,154721601,154721601), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(154853148,154853148,156440507,156440507), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(156841528,156841528,179857212,179857212), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(180691185,180691185,181000239,181000239), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(182105712,182105712,185834275,185834275), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet06.01 & pet06.02
chr_list <- c("P06", "N06.16", "C06.16", "C06")
y_list1 <- c(40,25,0,0)
y_list2 <- c(130,100,60,100)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A06", map2_chr == chr_list[i], map1_pos>15000000, map1_pos<90000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(15000000,90000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet06.01 & pet06.02")
  polygon(c(56219803,56219803,78749603,78749603), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(19377671,19377671,56046922,56046922), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#645084", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet06.03 
chr_list <- c("P06", "N06.16", "C06.16", "C06")
y_list1 <- c(20,0,0,0)
y_list2 <- c(60,40,10,20)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A06", map2_chr == chr_list[i], map1_pos>6000000, map1_pos<13000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(6000000,13000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet06.03")
  polygon(c(8592837,8592837,10582776,10582776), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet07.01 
chr_list <- c("P07.4", "N07", "C07", "C07")
y_list1 <- c(50,10,10,0)
y_list2 <- c(130,70,60,80)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A07", map2_chr == chr_list[i], map1_pos>100000000, map1_pos<150000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(100000000,150000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet07.01")
  polygon(c(116413080,116413080,130545045,130545045), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet08.01 
chr_list <- c("P08", "N08", "C08", "C08")
y_list1 <- c(70,50,60,30)
y_list2 <- c(140,130,110,85)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A08", map2_chr == chr_list[i], map1_pos>70000000, map1_pos<110000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(70000000,110000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet08.01")
  polygon(c(86786672,86786672,94581640,94581640), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet09.01
chr_list <- c("P09", "N09", "C09", "C09")
y_list1 <- c(0,0,0,10)
y_list2 <- c(80,160,100,110)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A09", map2_chr == chr_list[i], map1_pos>100000000, map1_pos<150000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(100000000,150000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet09.01")
  polygon(c(105305324,105305324,123481124,123481124), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(128805928,128805928,140852908,140852908), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet10.01
chr_list <- c("P10", "N10", "C10", "C10")
y_list1 <- c(0,0,0,0)
y_list2 <- c(28,19,18,12)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A10", map2_chr == chr_list[i], map1_pos>0, map1_pos<20000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,20000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet10.01")
  polygon(c(293669,293669,299675,299675), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(384092,384092,499852,499852), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(1235691,1235691,1268588,1268588), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(1716092,1716092,1776585,1776585), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(1845471,1845471,14398715,14398715), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet11.01
chr_list <- c("P11", "N11", "C11", "C11")
y_list1 <- c(0,0,0,0)
y_list2 <- c(70,80,20,25)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A11", map2_chr == chr_list[i], map1_pos>0, map1_pos<90000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,90000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet11.01")
  polygon(c(3320576,3320576,65747566,65747566), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet12.01
chr_list <- c("P17", "N17", "C17", "C17")
y_list1 <- c(55,60,58,55)
y_list2 <- c(80,140,72,90)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A12", map2_chr == chr_list[i], map1_pos>80000000, map1_pos<180000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(80000000,180000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet12.01")
  polygon(c(89671701,89671701,95073326,95073326), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(98049349,98049349,98177109,98177109), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(100451519,100451519,111696502,111696502), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(155860790,155860790,164385696,164385696), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet13.01
chr_list <- c("P13", "N13", "C13", "C13")
y_list1 <- c(70,85,40,40)
y_list2 <- c(110,110,100,80)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A13", map2_chr == chr_list[i], map1_pos>140000000, map1_pos<155000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(140000000,155000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet13.01")
  polygon(c(147947071,147947071,149109754,149109754), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet14.01 & pet14.2
chr_list <- c("P14", "N14", "C14", "C14")
y_list1 <- c(0,0,0,0)
y_list2 <- c(150,150,120,120)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A14", map2_chr == chr_list[i], map1_pos>70000000, map1_pos<200000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(70000000,200000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet14.01 & pet14.02")
  polygon(c(98475810,98475810,171068709,171068709), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(135376473,135376473,144002322,144002322), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#645084", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet16.01 & pet16.02
chr_list <- c("P12.16", "N06.16", "C06.16", "C06")
chr_list2 <- c("P16.17", "N16", "C16.17", "C16.17")
chr_list3 <- c("P17", "N17", "C17", "C17")
y_list1 <- c(60,0,20,50)
y_list2 <- c(170,180,100,130)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A16", map2_chr == chr_list[i], map1_pos>0, map1_pos<220000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,220000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet16.01 & pet16.02")
  temp2 <- pet_map_list[[i]] %>% filter(map1_chr == "A16", map2_chr == chr_list2[i], map1_pos>0, map1_pos<220000000)
  points(temp2$map2_pos ~ temp2$map1_pos, col = alpha(color_palette[temp2$color + col_list[i]], 0.6), pch = 17, cex = 1.5)
  temp3 <- pet_map_list[[i]] %>% filter(map1_chr == "A16", map2_chr == chr_list3[i], map1_pos>0, map1_pos<220000000)
  points(temp3$map2_pos ~ temp3$map1_pos, col = alpha(color_palette[temp3$color + col_list[i]], 0.6), pch = 18, cex = 1.5)
  polygon(c(6481680,6481680,10666924,10666924), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#645084", 0.2), border = NA)
  polygon(c(6566825,6566825,23262121,23262121), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(206419062,206419062,210181068,210181068), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(212077689,212077689,212145071,212145071), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(214469811,214469811,215512235,215512235), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(215796657,215796657,215824755,215824755), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(216448563,216448563,216569957,216569957), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(216654211,216654211,217103619,217103619), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(217608231,217608231,217617390,217617390), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(218176281,218176281,218224144,218224144), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(34379628,34379628,45760629,45760629), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(69141022,69141022,69300799,69300799), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(90030717,90030717,90927508,90927508), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(108631478,108631478,112509280,112509280), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(118370912,118370912,118376238,118376238), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(118617269,118617269,122506810,122506810), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(123753242,123753242,123760902,123760902), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(137513538,137513538,147855628,147855628), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), xlim = c(0,4), yaxt = 'n', xaxt='n', bty="n", pch = 15)
  points(temp2$map2_pos ~ rep(2, nrow(temp2)), col = alpha(color_palette[temp2$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), xlim = c(0,4), yaxt = 'n', xaxt='n', bty="n", pch = 15)
  points(temp3$map2_pos ~ rep(3, nrow(temp3)), col = alpha(color_palette[temp3$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), xlim = c(0,4), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet17.01
chr_list <- c("P17", "N17", "C17", "C17")
y_list1 <- c(0,0,0,0)
y_list2 <- c(45,80,60,40)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A17", map2_chr == chr_list[i], map1_pos>3000000, map1_pos<30000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(3000000,30000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet17.01")
  polygon(c(11688320,11688320,22678243,22678243), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), xlim = c(0,4), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet17.02
chr_list <- c("P17", "N17", "C17", "C17")
y_list1 <- c(20,50,30,0)
y_list2 <- c(80,100,70,80)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A17", map2_chr == chr_list[i], map1_pos>30000000, map1_pos<80000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(30000000,80000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet17.02")
  polygon(c(39484071,39484071,68196212,68196212), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), xlim = c(0,4), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# pet17.03 & pet17.04
chr_list <- c("P16.17", "N16", "C16.17", "C16.17")
y_list1 <- c(0,0,0,0)
y_list2 <- c(50,200,20,15)
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A17", map2_chr == chr_list[i], map1_pos>170000000, map1_pos<210000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(170000000,210000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "pet17.03 & pet17.04")
  polygon(c(188458401,188458401,203655689,203655689), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(203677179,203677179,203718262,203718262), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(203770587,203770587,203775207,203775207), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(204176880,204176880,204410173,204410173), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(204714679,204714679,204747520,204747520), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(194134669,194134669,203768215,203768215), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#645084", 0.2), border = NA)
  polygon(c(204279632,204279632,204409363,204409363), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#645084", 0.2), border = NA)
  polygon(c(204693877,204693877,204751769,204751769), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#645084", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + col_list[i]], 0.6),
       ylim = c(y_list1[i], y_list2[i]), xlim = c(0,4), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}


######################################
### --- Final argophyllus plots --- ###
######################################

layout(matrix(c(1,1,1,1,1,2,3,3,3,3,3,4,5,5,5,5,5,6), 3, 6, byrow = TRUE))
par(mai=c(0.2,0.2,0.2,0.2), mar=c(2, 4, 1.5, 1.5))
arg_map_list <- list(arg1_map, arg2_map, arg3_map)
# inversion, argophyllus
color_palette <- c("#953071", "#338732")

# arg03.01 & arg03.02
chr_list <- c("ARG03", "C03", "C03")
y_list1 <- c(0,0,0)
y_list2 <- c(50,50,50)
for (i in 1:3){
  temp <- arg_map_list[[i]] %>% filter(map1_chr == "A03", map2_chr == chr_list[i], map1_pos<150000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + 1], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,150000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "arg03.01 & arg03.02")
  polygon(c(42210710,42210710,42457211,42457211), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(43315730,43315730,44156163,44156163), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(44521626,44521626,78204102,78204102), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(78867692,78867692,80067728,80067728), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(80123927,80123927,83317909,83317909), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(87571939,87571939,89113195,89113195), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(98713561,98713561,98761502,98761502), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(28583609,28583609,42210629,42210629), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#645084", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + 1], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# arg06.01 & arg06.02
chr_list <- c("ARG06_15", "C06", "C06")
y_list1 <- c(0,0,0)
y_list2 <- c(100,100,100)
for (i in 1:3){
  temp <- arg_map_list[[i]] %>% filter(map1_chr == "A06", map2_chr == chr_list[i], map1_pos>50000000, map1_pos<205000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + 1], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(50000000,205000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "arg06.01 & arg06.02")
  polygon(c(130156472,130156472,153317127,153317127), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(153577387,153577387,154165841,154165841), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(154632762,154632762,154682381,154682381), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(154786258,154786258,155011244,155011244), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(155018155,155018155,155128857,155128857), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(125435746,125435746,130156455,130156455), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#645084", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + 1], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# arg10.01
chr_list <- c("ARG10", "C10", "C10")
y_list1 <- c(0,0,0)
y_list2 <- c(100,100,100)
for (i in 1:3){
  temp <- arg_map_list[[i]] %>% filter(map1_chr == "A10", map2_chr == chr_list[i], map1_pos>120000000, map1_pos<200000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + 1], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(120000000,200000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "arg10.01")
  polygon(c(161242043,161242043,190205715,190205715), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(190279711,190279711,190295665,190295665), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(190339415,190339415,190694486,190694486), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + 1], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# arg16.01
chr_list <- c("ARG16_12", "C12.16", "C12.16")
y_list1 <- c(0,0,0)
y_list2 <- c(50,50,50)
for (i in 1:3){
  temp <- arg_map_list[[i]] %>% filter(map1_chr == "A16", map2_chr == chr_list[i], map1_pos>190000000, map1_pos<220000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + 1], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(190000000,220000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "arg16.01")
  polygon(c(202452179,202452179,204242807,204242807), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(204537714,204537714,204629966,204629966), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(205197484,205197484,205209334,205209334), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + 1], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}


##################################
### --- Final annuus plots --- ###
##################################

layout(matrix(c(1,1,1,1,1,2,3,3,3,3,3,4,5,5,5,5,5,6), 3, 6, byrow = TRUE))
par(mai=c(0.2,0.2,0.2,0.2), mar=c(2, 4, 1.5, 1.5))
ann_map_list <- list(ann1_map, ann2_map, ann3_map)
# inversion, annuus
color_palette <- c("#953071", "#FFB31C")

# ann01.01
y_list1 <- c(0,0,0)
y_list2 <- c(20,20,20)
for (i in 1:3){
  temp <- ann_map_list[[i]] %>% filter(map1_chr == "A01", map2_chr == 1, map1_pos>0, map1_pos<25000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + 1], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,25000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "ann01.01")
  polygon(c(6199548,6199548,9699260,9699260), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + 1], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# ann05.01
y_list1 <- c(20,20,20)
y_list2 <- c(100,100,100)
for (i in 1:3){
  temp <- ann_map_list[[i]] %>% filter(map1_chr == "A05", map2_chr == 5, map1_pos>120000000, map1_pos<200000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + 1], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(120000000,200000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "ann05.01")
  polygon(c(148181142,148181142,177442806,177442806), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + 1], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# ann11.01
y_list1 <- c(0,0,0)
y_list2 <- c(80,80,80)
for (i in 1:3){
  temp <- ann_map_list[[i]] %>% filter(map1_chr == "A11", map2_chr == 11, map1_pos>0, map1_pos<80000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + 1], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,80000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "ann11.01")
  polygon(c(19370443,19370443,53975915,53975915), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + 1], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# ann13.01 & ann13.02
y_list1 <- c(0,0,0)
y_list2 <- c(100,100,100)
for (i in 1:3){
  temp <- ann_map_list[[i]] %>% filter(map1_chr == "A13", map2_chr == 13, map1_pos>0, map1_pos<180000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + 1], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,180000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "ann13.01 & ann13.02")
  polygon(c(8646123,8646123,109475871,109475871), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(139223259,139223259,157756098,157756098), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#645084", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + 1], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# ann14.01 & ann14.02
y_list1 <- c(0,0,0)
y_list2 <- c(60,60,60)
for (i in 1:3){
  temp <- ann_map_list[[i]] %>% filter(map1_chr == "A14", map2_chr == 14, map1_pos>70000000, map1_pos<160000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + 1], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(70000000,160000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "ann14.01 & ann14.02")
  polygon(c(101455802,101455802,129356206,129356206), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  polygon(c(129356233,129356233,133235392,133235392), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#645084", 0.2), border = NA)
  polygon(c(133459531,133459531,134701278,134701278), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#645084", 0.2), border = NA)
  polygon(c(135046321,135046321,135377910,135377910), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#645084", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + 1], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# ann15.01
y_list1 <- c(20,20,20)
y_list2 <- c(100,100,100)
for (i in 1:3){
  temp <- ann_map_list[[i]] %>% filter(map1_chr == "A15", map2_chr == 15, map1_pos>50000000, map1_pos<190000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + 1], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(50000000,190000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "ann15.01")
  polygon(c(103772702,103772702,176634307,176634307), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + 1], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# ann16.01
y_list1 <- c(0,0,0)
y_list2 <- c(80,80,80)
for (i in 1:3){
  temp <- ann_map_list[[i]] %>% filter(map1_chr == "A16", map2_chr == 16, map1_pos>0, map1_pos<50000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + 1], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,50000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "ann16.01")
  polygon(c(10542720,10542720,23486443,23486443), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + 1], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# ann16.02
y_list1 <- c(30,30,30)
y_list2 <- c(80,80,80)
for (i in 1:3){
  temp <- ann_map_list[[i]] %>% filter(map1_chr == "A16", map2_chr == 16, map1_pos>120000000, map1_pos<180000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + 1], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(120000000,180000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "ann16.02")
  polygon(c(147445328,147445328,159006308,159006308), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + 1], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

# ann17.01
y_list1 <- c(50,50,50)
y_list2 <- c(110,110,110)
for (i in 1:3){
  temp <- ann_map_list[[i]] %>% filter(map1_chr == "A17", map2_chr == 17, map1_pos>180000000, map1_pos<210000000)
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + 1], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(180000000,210000000), ylim = c(y_list1[i], y_list2[i]), las = 1, cex = 1.5, main = "ann17.01")
  polygon(c(188468191,188468191,197519081,197519081), c(y_list1[i],y_list2[i],y_list2[i],y_list1[i]), col=alpha("#953071", 0.2), border = NA)
  plot(temp$map2_pos ~ rep(1, nrow(temp)), col = alpha(color_palette[temp$color + 1], 0.6),
       ylim = c(y_list1[i], y_list2[i]), yaxt = 'n', xaxt='n', bty="n", pch = 15)
}

############################
### --- Old pet maps --- ###
############################

layout(matrix(1:4, 4, 1, byrow = FALSE))
par(mai=c(0.2,0.2,0.2,0.2), mar=c(2, 4, 1.5, 1.5))
pet_map_list <- list(pet1_map, pet2_map, pet3_map, pet4_map)
# inversion, fallax, inversion, petiolaris
color_palette <- c("#645084", "#447499", "#645084", "#4BB7B7")
col_list <- c(3,1,1,1)

# pet01.01
chr_list <- c("P01", "N01", "C01", "C01")
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A01", map2_chr == chr_list[i])
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,160000000), ylim = NULL, las = 1, cex = 1.5, main = "pet01.01")
}

# pet05.01
chr_list <- c("P05", "N05", "C05", "C05")
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A05", map2_chr == chr_list[i])
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,182000000), ylim = NULL, las = 1, cex = 1.5, main = "pet05.01")
}

# pet06.01 & pet06.02 & pet06.03
chr_list <- c("P06", "N06.16", "C06.16", "C06")
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A06", map2_chr == chr_list[i])
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,150000000), ylim = NULL, las = 1, cex = 1.5, main = "pet06.01 & pet06.02 & pet06.03")
  abline(v=c(56219803, 78749603), lty=1)
  abline(v=c(19377671, 56046922), lty=2)
  abline(v=c(8592837, 10582776), lty=3)
}

# pet07.01
chr_list <- c("P07.4", "N07", "C07", "C07")
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A07", map2_chr == chr_list[i])
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,170000000), ylim = NULL, las = 1, cex = 1.5, main = "pet07.01")
}

# pet08.01
chr_list <- c("P08", "N08", "C08", "C08")
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A08", map2_chr == chr_list[i])
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,175000000), ylim = NULL, las = 1, cex = 1.5, main = "pet08.01")
}

# pet09.01
chr_list <- c("P09", "N09", "C09", "C09")
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A09", map2_chr == chr_list[i])
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,200000000), ylim = NULL, las = 1, cex = 1.5, main = "pet09.01")
}

# pet10.01
chr_list <- c("P10", "N10", "C10", "C10")
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A10", map2_chr == chr_list[i])
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,190000000), ylim = NULL, las = 1, cex = 1.5, main = "pet10.01")
}

# pet11.01
chr_list <- c("P11", "N11", "C11", "C11")
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A11", map2_chr == chr_list[i])
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,200000000), ylim = NULL, las = 1, cex = 1.5, main = "pet11.01")
}

# pet12.01
chr_list <- c("P17", "N17", "C17", "C17")
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A12", map2_chr == chr_list[i])
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(50000000,180000000), ylim = NULL, las = 1, cex = 1.5, main = "pet12.01")
}

# pet13.01
chr_list <- c("P13", "N13", "C13", "C13")
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A13", map2_chr == chr_list[i])
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,200000000), ylim = NULL, las = 1, cex = 1.5, main = "pet13.01")
}

# pet14.01
chr_list <- c("P14", "N14", "C14", "C14")
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A14", map2_chr == chr_list[i])
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,200000000), ylim = NULL, las = 1, cex = 1.5, main = "pet14.01 & pet14.02")
  abline(v=c(135376473, 144002322))
}

# pet16.01 & pet16.02
chr_list <- c("P12.16", "N06.16", "C06.16", "C06")
chr_list2 <- c("P16.17", "N16", "C16.17", "C16.17")
chr_list3 <- c("P17", "N17", "C17", "C17")
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A16", map2_chr == chr_list[i])
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,220000000), ylim = c(0,145), las = 1, cex = 1.5, main = "pet16.01 & pet16.02")
  temp2 <- pet_map_list[[i]] %>% filter(map1_chr == "A16", map2_chr == chr_list2[i])
  points(temp2$map2_pos ~ temp2$map1_pos, col = alpha(color_palette[temp2$color + col_list[i]], 0.6), pch = 17, cex = 1.5)
  temp3 <- pet_map_list[[i]] %>% filter(map1_chr == "A16", map2_chr == chr_list3[i])
  points(temp3$map2_pos ~ temp3$map1_pos, col = alpha(color_palette[temp3$color + col_list[i]], 0.6), pch = 18, cex = 1.5)
  abline(v=c(6481680, 10666924))
}

# pet17.01-pet17.4
chr_list <- c("P17", "N17", "C17", "C17")
chr_list2 <- c("P16.17", "N16", "C16.17", "C16.17")
for (i in 1:4){
  temp <- pet_map_list[[i]] %>% filter(map1_chr == "A17", map2_chr == chr_list[i])
  plot(temp$map2_pos ~ temp$map1_pos, col = alpha(color_palette[temp$color + col_list[i]], 0.6), 
       xlab = NULL, ylab = "cM", pch = 16, xlim = c(0,210000000), ylim = c(0,145), las = 1, cex = 1.5, main = "pet17.01 - pet17.04")
  temp2 <- pet_map_list[[i]] %>% filter(map1_chr == "A17", map2_chr == chr_list2[i])
  points(temp2$map2_pos ~ temp2$map1_pos, col = alpha(color_palette[temp2$color + col_list[i]], 0.6), pch = 17, cex = 1.5)
  abline(v=c(11688320, 22678243), lty=1)
  abline(v=c(39484071, 68196212), lty=2)
  abline(v=c(194134669, 203768215, 204279632, 204409363, 204693877, 204751769), lty=3)
}
