library(tidyverse)

data <- tibble(pos1=numeric(),pos2=numeric(),link=numeric())
circle_size <- 4
circle_y <- -0.5

for (i in 1:30){
  for (j in 1:30){
    distance <- abs(j - i)
    link <- 5 - distance
    if (link < 0){
      link = 0
    }
    tmp <- tibble(pos1=i,pos2=j,link=link)
    data <- rbind(data,tmp)
  }
}
wild_type_data <- data %>% 
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
  rename(original_link = link) 


pdf("HiC/inversion_examples.pdf",height=3,width=4)

data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(color_value = pos1) %>%
  select(pos1,color_value) %>%
  unique() -> color_dots

data %>% 
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
  ggplot(.,aes()) + geom_polygon(aes(x=x,y=y,fill=link,group = bin_id),size=0) +
  geom_point(data=color_dots, aes(x=pos1,y=circle_y,color=color_value),size=circle_size,show.legend = F) +
  scale_colour_distiller(palette = "Greens",) +
  scale_fill_gradient(low = "#377EB8",
                       high = "#E41A1C", name="HiC interactions") +

  theme_linedraw() + ylab("") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        panel.border = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

####One small inversion

data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(color_value = pos1) %>%
  mutate(pos1 = case_when(pos1 >=5 & pos1 <= 10 ~ (10 - pos1 + 5),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=5 & pos2 <= 10 ~ (10 - pos2+ 5),
                          TRUE ~ pos2)) %>% 
  select(pos1,color_value) %>%
  unique() -> color_dots

data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(pos1 = case_when(pos1 >=5 & pos1 <= 10 ~ (10 - pos1 + 5),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=5 & pos2 <= 10 ~ (10 - pos2+ 5),
                          TRUE ~ pos2)) %>% 
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
  ggplot(.,aes()) + geom_polygon(aes(x=x,y=y,fill=link,group = bin_id),size=0) +
  geom_point(data=color_dots, aes(x=pos1,y=circle_y,color=color_value),size=circle_size,show.legend = F) +
  scale_colour_distiller(palette = "Greens",) +
  scale_fill_gradient(low = "#377EB8",
                      high = "#E41A1C", name="HiC interactions") +
  
  theme_linedraw() + ylab("") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        panel.border = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(pos1 = case_when(pos1 >=5 & pos1 <= 10 ~ (10 - pos1 + 5),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=5 & pos2 <= 10 ~ (10 - pos2+ 5),
                          TRUE ~ pos2)) %>% 
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
  inner_join(.,wild_type_data) %>%
  mutate(interaction_dif = link - original_link ) %>%
  ggplot(.,aes()) + geom_polygon(aes(x=x,y=y,fill=interaction_dif,group = bin_id),size=0) +
  geom_point(data=color_dots, aes(x=pos1,y=circle_y,color=color_value),size=circle_size,show.legend = F) +
  scale_colour_distiller(palette = "Greens") +
  scale_fill_gradient2(low = "#377EB8", mid = "#f2f2f2",
                       high = "#E41A1C", midpoint = 0,name="HiC interaction\ndifference",limits=c(-4,4)) +
  
  theme_linedraw() + ylab("") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        panel.border = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

#####One big inversions
data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(color_value = pos1) %>%
  mutate(pos1 = case_when(pos1 >=9 & pos1 <= 25 ~ (25 - pos1 + 9),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=9 & pos2 <= 25 ~ (25 - pos2+ 9),
                          TRUE ~ pos2)) %>% 
  select(pos1,color_value) %>%
  unique() -> color_dots


data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(pos1 = case_when(pos1 >=9 & pos1 <= 25 ~ (25 - pos1 + 9),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=9 & pos2 <= 25 ~ (25 - pos2+ 9),
                          TRUE ~ pos2)) %>% 
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
  ggplot(.,aes()) + geom_polygon(aes(x=x,y=y,fill=link,group = bin_id),size=0) +
  geom_point(data=color_dots, aes(x=pos1,y=circle_y,color=color_value),size=circle_size,show.legend = F) +
  scale_colour_distiller(palette = "Greens") +
  scale_fill_gradient(low = "#377EB8",
                      high = "#E41A1C", name="HiC interactions") +
  
  theme_linedraw() + ylab("") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        panel.border = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 


data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(pos1 = case_when(pos1 >=9 & pos1 <= 25 ~ (25 - pos1 + 9),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=9 & pos2 <= 25 ~ (25 - pos2+ 9),
                          TRUE ~ pos2)) %>% 
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
inner_join(.,wild_type_data) %>%
  mutate(interaction_dif = link - original_link ) %>%
  ggplot(.,aes()) + geom_polygon(aes(x=x,y=y,fill=interaction_dif,group = bin_id),size=0) +
  geom_point(data=color_dots, aes(x=pos1,y=circle_y,color=color_value),size=circle_size,show.legend = F) +
  scale_colour_distiller(palette = "Greens") +
  scale_fill_gradient2(low = "#377EB8", mid = "#f2f2f2",
                       high = "#E41A1C", midpoint = 0,name="HiC interaction\ndifference",limits=c(-4,4)) +
  
  theme_linedraw() + ylab("") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        panel.border = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 




#####Two inversions

data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(color_value = pos1) %>%
  mutate(pos1 = case_when(pos1 >=5 & pos1 <= 10 ~ (10 - pos1 + 5),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=5 & pos2 <= 10 ~ (10 - pos2+ 5),
                          TRUE ~ pos2)) %>% 
  mutate(pos1 = case_when(pos1 >=9 & pos1 <= 25 ~ (25 - pos1 + 9),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=9 & pos2 <= 25 ~ (25 - pos2+ 9),
                          TRUE ~ pos2)) %>%
  select(pos1,color_value) %>%
  unique() -> color_dots


data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(pos1 = case_when(pos1 >=5 & pos1 <= 10 ~ (10 - pos1 + 5),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=5 & pos2 <= 10 ~ (10 - pos2+ 5),
                          TRUE ~ pos2)) %>% 
  mutate(pos1 = case_when(pos1 >=9 & pos1 <= 25 ~ (25 - pos1 + 9),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=9 & pos2 <= 25 ~ (25 - pos2+ 9),
                          TRUE ~ pos2)) %>% 
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
  ggplot(.,aes()) + geom_polygon(aes(x=x,y=y,fill=link,group = bin_id),size=0) +
  geom_point(data=color_dots, aes(x=pos1,y=circle_y,color=color_value),size=circle_size,show.legend = F) +
  scale_colour_distiller(palette = "Greens") +
  scale_fill_gradient(low = "#377EB8",
                      high = "#E41A1C", name="HiC interactions") +
  
  theme_linedraw() + ylab("") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        panel.border = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 




data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(pos1 = case_when(pos1 >=5 & pos1 <= 10 ~ (10 - pos1 + 5),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=5 & pos2 <= 10 ~ (10 - pos2+ 5),
                          TRUE ~ pos2)) %>% 
  mutate(pos1 = case_when(pos1 >=9 & pos1 <= 25 ~ (25 - pos1 + 9),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=9 & pos2 <= 25 ~ (25 - pos2+ 9),
                          TRUE ~ pos2)) %>% 
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
  inner_join(.,wild_type_data) %>%
  mutate(interaction_dif = link - original_link ) %>%
  ggplot(.,aes()) + geom_polygon(aes(x=x,y=y,fill=interaction_dif,group = bin_id),size=0) +
  geom_point(data=color_dots, aes(x=pos1,y=circle_y,color=color_value),size=circle_size,show.legend = F) +
  scale_colour_distiller(palette = "Greens") +
  scale_fill_gradient2(low = "#377EB8", mid = "#f2f2f2",
                       high = "#E41A1C", midpoint = 0,name="HiC interaction\ndifference",limits=c(-4,4)) +
  
  theme_linedraw() + ylab("") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        panel.border = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

dev.off()
##############################
####For pet09.01
##############################
pdf("HiC/inversion_examples_forpaper.pdf",height=3,width=6)


data <- tibble(pos1=numeric(),pos2=numeric(),link=numeric())
circle_size <- 4
circle_y <- -0.5

for (i in 90:160){
  for (j in 90:160){
    distance <- abs(j - i)
    link <- 5 - distance
    if (link < 0){
      link = 0
    }
    tmp <- tibble(pos1=i,pos2=j,link=link)
    data <- rbind(data,tmp)
  }
}
wild_type_data <- data %>% 
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(pos1 = case_when(pos1 >=124 & pos1 <= 140 ~ (124 - pos1 + 140),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=124 & pos2 <= 140 ~ (124 - pos2 + 140),
                          TRUE ~ pos2)) %>% 
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
  rename(original_link = link) 


data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(color_value = pos1) %>%
  mutate(pos1 = case_when(pos1 >=105 & pos1 <= 135 ~ (105 - pos1 + 135),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=105 & pos2 <= 135 ~ (105 - pos2 + 135),
                          TRUE ~ pos2)) %>% 
  mutate(pos1 = case_when(pos1 >=124 & pos1 <= 140 ~ (124 - pos1 + 140),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=124 & pos2 <= 140 ~ (124 - pos2 + 140),
                          TRUE ~ pos2)) %>% 
  select(pos1,color_value) %>%
  unique() -> color_dots



data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(pos1 = case_when(pos1 >=105 & pos1 <= 135 ~ (105 - pos1 + 135),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=105 & pos2 <= 135 ~ (105 - pos2 + 135),
                          TRUE ~ pos2)) %>% 
  mutate(pos1 = case_when(pos1 >=124 & pos1 <= 140 ~ (124 - pos1 + 140),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=124 & pos2 <= 140 ~ (124 - pos2 + 140),
                          TRUE ~ pos2)) %>% 
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
  inner_join(.,wild_type_data) %>%
  mutate(interaction_dif = link - original_link ) %>%
  ggplot(.,aes()) + geom_polygon(aes(x=x,y=y,fill=interaction_dif,group = bin_id),size=0) +
  geom_point(data=color_dots, aes(x=pos1,y=circle_y,color=color_value),size=circle_size,show.legend = F) +
  scale_colour_distiller(palette = "Greens") +
  scale_fill_gradient2(low = "#377EB8", mid = "#f2f2f2",
                       high = "#E41A1C", midpoint = 0,name="HiC interaction\ndifference",limits=c(-4,4)) +
  
  theme_linedraw() + ylab("") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        panel.border = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  annotate("segment",x=105,y=0,
           yend=17.5,xend=122.5,
           alpha=1,color="black",size=0.1) +
  annotate("segment",x=140,y=0,
           yend=17.5,xend=122.5,
           alpha=1,color="black",size=0.1) +
  annotate("segment",x=105,y=0,
         yend=9.5,xend=114.5,
         alpha=1,color="black",size=0.1,linetype="dotted") +
  annotate("segment",x=124,y=0,
           yend=9.5,xend=114.5,
           alpha=1,color="black",size=0.1,linetype="dotted") +
  annotate("segment",x=128,y=0,
           yend=6,xend=134,
           alpha=1,color="black",size=0.1,linetype="dotted") +
  annotate("segment",x=140,y=0,
           yend=6,xend=134,
           alpha=1,color="black",size=0.1,linetype="dotted") +
  ggtitle("pet09.01")


####Pet05.01


data <- tibble(pos1=numeric(),pos2=numeric(),link=numeric())
circle_size <- 4
circle_y <- -0.5

for (i in 140:187){
  for (j in 140:187){
    distance <- abs(j - i)
    link <- 5 - distance
    if (link < 0){
      link = 0
    }
    tmp <- tibble(pos1=i,pos2=j,link=link)
    data <- rbind(data,tmp)
  }
}

wild_type_data <- data %>% 
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(pos1 = case_when(pos1 >=181 & pos1 <= 186 ~ (181 - pos1 + 186),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=181 & pos2 <= 186 ~ (181 - pos2 + 186),
                          TRUE ~ pos2)) %>% 
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
  rename(original_link = link) 



data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(color_value = pos1) %>%
  mutate(pos1 = case_when(pos1 >=157 & pos1 <= 184 ~ (157 - pos1 + 184),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=157 & pos2 <= 184 ~ (157 - pos2 + 184),
                          TRUE ~ pos2)) %>% 
  mutate(pos1 = case_when(pos1 >=181 & pos1 <= 185 ~ (181 - pos1 + 185),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=181 & pos2 <= 185 ~ (181 - pos2 + 185),
                          TRUE ~ pos2)) %>% 
  select(pos1,color_value) %>%
  unique() -> color_dots



data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(pos1 = case_when(pos1 >=157 & pos1 <= 184 ~ (157 - pos1 + 184),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=157 & pos2 <= 184 ~ (157 - pos2 + 184),
                          TRUE ~ pos2)) %>% 
  mutate(pos1 = case_when(pos1 >=181 & pos1 <= 185 ~ (181 - pos1 + 185),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=181 & pos2 <= 185 ~ (181 - pos2 + 185),
                          TRUE ~ pos2)) %>% 
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
  inner_join(.,wild_type_data) %>%
  mutate(interaction_dif = link - original_link ) %>%
  ggplot(.,aes()) + geom_polygon(aes(x=x,y=y,fill=interaction_dif,group = bin_id),size=0) +
  geom_point(data=color_dots, aes(x=pos1,y=circle_y,color=color_value),size=circle_size,show.legend = F) +
  scale_colour_distiller(palette = "Greens") +
  scale_fill_gradient2(low = "#377EB8", mid = "#f2f2f2",
                       high = "#E41A1C", midpoint = 0,name="HiC interaction\ndifference",limits=c(-4,4)) +
  
  theme_linedraw() + ylab("") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        panel.border = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  annotate("segment",x=154,y=0,
           yend=15.5,xend=169.5,
           alpha=1,color="black",size=0.1) +
  annotate("segment",x=186,y=0,
           yend=15.5,xend=169.5,
           alpha=1,color="black",size=0.1) +
  annotate("segment",x=154,y=0,
           yend=13.5,xend=167.5,
           alpha=1,color="black",size=0.1,linetype="dotted") +
  annotate("segment",x=181,y=0,
           yend=13.5,xend=167.5,
           alpha=1,color="black",size=0.1,linetype="dotted") +
  annotate("segment",x=182,y=0,
           yend=2,xend=184,
           alpha=1,color="black",size=0.1,linetype="dotted") +
  annotate("segment",x=186,y=0,
           yend=2,xend=184,
           alpha=1,color="black",size=0.1,linetype="dotted") +
  ggtitle("pet05.01")


#####Pet17.03



data <- tibble(pos1=numeric(),pos2=numeric(),link=numeric())
circle_size <- 4
circle_y <- -0.5

for (i in 169:205){
  for (j in 169:205){
    distance <- abs(j - i)
    link <- 5 - distance
    if (link < 0){
      link = 0
    }
    tmp <- tibble(pos1=i,pos2=j,link=link)
    data <- rbind(data,tmp)
  }
}
wild_type_data <- data %>% 
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
  rename(original_link = link) 


data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(color_value = pos1) %>%
  mutate(pos1 = case_when(pos1 >=188 & pos1 <= 204 ~ (188 - pos1 + 204),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=188 & pos2 <= 204 ~ (188 - pos2 + 204),
                          TRUE ~ pos2)) %>% 
  select(pos1,color_value) %>%
  unique() -> color_dots



data %>%
  mutate(pos1 = as.numeric(pos1),
         pos2 = as.numeric(pos2)) %>%
  mutate(pos1 = case_when(pos1 >=188 & pos1 <= 204 ~ (188 - pos1 + 204),
                          TRUE ~ pos1),
         pos2 = case_when(pos2 >=188 & pos2 <= 204 ~ (188 - pos2 + 204),
                          TRUE ~ pos2)) %>% 
  filter(pos2 >= pos1) %>%
  mutate(x = as.numeric((pos1 + pos2)/2),
         y=as.numeric((pos2-pos1)/2)) %>%
  group_by(y,x,link) %>%
  expand(count = seq(1:4)) %>%
  mutate(bin_id = paste0(x,"_",y)) %>%
  ungroup() %>%
  mutate(y = case_when(count == 2 ~ y + 0.5,
                       count == 4 ~ y - 0.5,
                       TRUE ~ y),
         x = case_when(count == 1 ~ x - 0.5,
                       count == 3 ~ x + 0.5,
                       TRUE ~ x)) %>%
  mutate(y = case_when(y < 0 ~ 0,
                       TRUE ~ y)) %>%
  inner_join(.,wild_type_data) %>%
  mutate(interaction_dif = link - original_link ) %>%
  ggplot(.,aes()) + geom_polygon(aes(x=x,y=y,fill=interaction_dif,group = bin_id),size=0) +
  geom_point(data=color_dots, aes(x=pos1,y=circle_y,color=color_value),size=circle_size,show.legend = F) +
  scale_colour_distiller(palette = "Greens") +
  scale_fill_gradient2(low = "#377EB8", mid = "#f2f2f2",
                       high = "#E41A1C", midpoint = 0,name="HiC interaction\ndifference",limits=c(-4,4)) +
  
  theme_linedraw() + ylab("") + xlab("Mbp") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        panel.border = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  annotate("segment",x=188,y=0,
           yend=8,xend=196,
           alpha=1,color="black",size=0.1) +
  annotate("segment",x=204,y=0,
           yend=8,xend=196,
           alpha=1,color="black",size=0.1) +
  ggtitle("pet17.03")


dev.off()

