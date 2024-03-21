# Written in 2023 by Lars Koch
# Contact: https://www.linkedin.com/in/lars-koch-93844023a/
# Publication: "ADAMTS12 promotes fibrosis by restructuring ECM to enable activation 
# and migration of injury-responsive fibroblasts"

library(ggplot2)
library(xlsx)
library(stringr)
library(dplyr)
library(mgsub)
library(ggpubr)
setwd("xxx") # Set your working directory here

# Load spot dataset retrieved from ImageJ Plugin TrackMate 
# https://imagej.net/plugins/trackmate/
spots <- rbind(readxl::read_excel("Input/Matrigel lci5 spots 277.xlsx"),
               readxl::read_excel("Input/Matrigel lci5 spots NG.xlsx"))

# Load track dataset from TrackMate
spot_track <- readxl::read_excel("Input/Matrigel lci5 tracks.xlsx")

# Arrange spot dataset
spots$id <- paste(str_split_fixed(spots$Name, pattern = " spot", n = 2)[,1],
                  spots$TRACK_ID)
spots <- spots[,c("id","POSITION_X", "POSITION_Y", "POSITION_T")]
spots$POSITION_X <- as.numeric(spots$POSITION_X)
spots$POSITION_Y <- as.numeric(spots$POSITION_Y)
spots$POSITION_T <- as.numeric((spots$POSITION_T))
spots <- na.omit(spots)

# Arrange track dataset
spot_track$id <- paste(str_split_fixed(spot_track$Name, pattern = " track", n = 2)[,1],
                   spot_track$TRACK_ID)
spot_track <- spot_track[, c("id", "TRACK_MEAN_SPEED", "NUMBER_SPOTS")]
spot_track$TRACK_MEAN_SPEED <- as.numeric(spot_track$TRACK_MEAN_SPEED)
spot_track$NUMBER_SPOTS <- as.numeric(spot_track$NUMBER_SPOTS)
spot_track <- na.omit(spot_track)

# Merge spot and track dataset to combine Average Speed and Position in one dataframe
spot_track <- merge(spots, spot_track, by = "id")

# Name groups
spot_track$group <- paste(str_split_fixed(spot_track$id, " ", n = 3)[,1], 
                          str_split_fixed(spot_track$id, " ", n = 3)[,2])
spot_track$group <- 
  mgsub(spot_track$group,
        c("NG CT", "NG TGFb", "277 CT", "277 TGFb"), 
        c("WT CT", "WT TGFb", "KO CT", "KO TGFb"))

# Arrange merged dataframe
spot_track$stim <- str_split_fixed(spot_track$group, " ",2)[,2]
spot_track$genotype <- str_split_fixed(spot_track$group, " ",2)[,1]
spot_track$genotype <- gsub("KO", "ADAMTS12-KO", x = spot_track$genotype)
spot_track$genotype <- factor(spot_track$genotype, levels = c("WT", "ADAMTS12-KO"))

# Convert speed from pixel/10min to um/h and position from pixel to um
spot_track$TRACK_MEAN_SPEED <- spot_track$TRACK_MEAN_SPEED*14.781605114
spot_track$POSITION_X <- spot_track$POSITION_X*2.46360085
spot_track$POSITION_Y <- spot_track$POSITION_Y*2.46360085

# Extract starting point for each track
start_point <- spot_track %>% group_by(id) %>% slice(which.min(POSITION_T))
start_point <- start_point[c("id", "POSITION_X", "POSITION_Y")]
colnames(start_point) <- c("id", "start_x", "start_y")

# Substracting starting point from each coordinate to center each track starting from 0,0
spot_track <- merge(spot_track, start_point)
spot_track$zero_x <- spot_track$POSITION_X - spot_track$start_x
spot_track$zero_y <- spot_track$POSITION_Y - spot_track$start_y

# Order by timepoint and group
spot_track <- spot_track[order(spot_track$id, spot_track$POSITION_T),]
spot_track <- transform(spot_track,
                        group=factor(group,levels=c("WT CT","WT TGFb","KO CT", "KO TGFb")))

########################## Fig. 4E - Trajectory Maps ###########################

ggplot(spot_track, aes(x = spot_track$zero_x, y = spot_track$zero_y, 
                       group = spot_track$id, color = spot_track$TRACK_MEAN_SPEED))+
  geom_hline(aes(yintercept = 0), size = 0.3)+
  geom_vline(aes(xintercept = 0), size = 0.3)+
  geom_path(size = 0.3)+    
  xlim(-450, 450)+
  ylim(-450, 450)+
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 20,
    limits = c(0,40))+
  xlab("x-axis/µm")+ 
  ylab("y-axis/µm")+
  labs(color = "speed\nµm/h")+
  theme_void()+
  theme(plot.title = element_blank(),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        axis.title.y = element_text(size =8, face= "plain", family = "Arial", angle = 90),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", angle = 0),
        legend.text = element_text(size =6, face= "plain", family = "Arial"),
        legend.title = element_text(size =8, face= "plain", family = "Arial"),
        legend.key.height = unit(0.1, "inches"),
        legend.key.width = unit(0.05, "inches"),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_line(colour = "grey", size = 0.05))+
  facet_grid(genotype ~ stim)+
  theme(strip.text.y = element_blank())

ggsave(paste0("Output/Fig. 4E - Trajectory Maps.tiff"), width = 2.5, height = 2.1, units = "in", dpi = 1000)

############################ Fig. 4E - Analysis ################################

# Weight speed of each track by th captured spots in the track
spot_track$spotsxspeed <- spot_track$NUMBER_SPOTS*spot_track$TRACK_MEAN_SPEED

# Summarize
spot_track$pic <- sub("^(\\w+\\s+[A-Za-z]+\\s+\\d+_\\d+).*", "\\1", spot_track$id)
spot_track <- group_by(.data = spot_track, pic)
Erg <- summarize(spot_track, sum(NUMBER_SPOTS), sum(spotsxspeed))
Erg$average_speed_per_frame <- Erg$`sum(spotsxspeed)`/Erg$`sum(NUMBER_SPOTS)`

writexl::write_xlsx(path = "Output/Fig. 4E - Analysis.xlsx", x = Erg)

########################### Supp Video 2 - Paths ###############################

spot_track$pow <- paste(str_split_fixed(spot_track$id, " ", 4)[,1], 
                        str_split_fixed(spot_track$id, " ", 4)[,2],
                        str_split_fixed(spot_track$id, " ", 4)[,3])

# Define ROIs to use in video
pics_for_vid <- c("277 CT 2_2", "277 TGFb 3_2", "NG CT 1_3", "NG TGFb 1_1")

for(i in seq(1,length(pics_for_vid),1)){
  # Subset only ROI
  vid <- spot_track[spot_track$pow == pics_for_vid[i],]
  # Plot
  ggplot(vid, aes(x = POSITION_X, y = -POSITION_Y, group = id, color = TRACK_MEAN_SPEED))+
    geom_path(size = 0.5)+
    scale_colour_gradient2(
      low = "darkblue",
      mid = "whitesmoke",
      high = "indianred",
      midpoint = 20,
      limits = c(0,40))+
    scale_x_continuous(expand = c(0, 0))+ 
    scale_y_continuous(expand = c(0, 0))+
    labs(color='Speed [µm/h]')+
    theme_void()+
    theme(legend.position = "none")
  
  ggsave(paste0("Output/Supp Video 2 - Paths/" , pics_for_vid[i], ".svg"), width = 4, height = 4)
}



