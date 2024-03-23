# Written in 2023 by Lars Koch
# Contact: https://www.linkedin.com/in/lars-koch-93844023a/
# Publication: "ADAMTS12 promotes fibrosis by restructuring ECM to enable activation 
# and migration of injury-responsive fibroblasts"

library(stringr)
library(dplyr)
library(ggplot2)
library(gganimate)
library(viridis) 
library(mgsub)
setwd("xxx") # Set your working directory here

# Import and arrange track datafrmae
dir <- "Input/CSV/TRACKS"
files <- paste(1:92, "t.csv")
df_list <- lapply(paste0(dir, "/", files), function(x) read.csv(x))

# Convert list into dataframe, with filename as ID
df_tracks <- data.frame()
for(i in seq(1, length(df_list), 1)){
  df_tracks <- rbind(df_tracks, data.frame(df_list[[i]], id = i))
}

df_tracks <- df_tracks[, c("id", "TRACK_INDEX", "TRACK_MEAN_SPEED", "NUMBER_SPOTS")]
df_tracks$TRACK_MEAN_SPEED <- as.numeric(df_tracks$TRACK_MEAN_SPEED)
df_tracks$NUMBER_SPOTS <- as.numeric(df_tracks$NUMBER_SPOTS)
df_tracks$NUMBER_EDGES <- df_tracks$NUMBER_SPOTS-1
df_tracks <- na.omit(df_tracks)

# Assign group
df_tracks$Genotype <- ifelse(df_tracks$id <= 32, "B2",
                             ifelse(df_tracks$id <= 64, "B3",
                                    ifelse(df_tracks$id <= 92, "B4", NA)))

df_tracks$Stim <- ifelse(df_tracks$id <= 16, "CT", 
                         ifelse(df_tracks$id <= 32, "TGFb",
                                ifelse(df_tracks$id <= 48, "CT",
                                       ifelse(df_tracks$id <= 64, "TGFb",
                                              ifelse(df_tracks$id <= 80, "CT",
                                                     ifelse(df_tracks$id <= 92, "TGFb", NA))))))

# Import and arrange spot datafrmae
dir <- "Input/CSV/SPOTS"
files <- paste(1:92, "S.csv")
df_list <- lapply(paste0(dir, "/", files), function(x) read.csv(x))

# Convert list into dataframe, with filename as ID
df_spots <- data.frame()
for(i in seq(1, length(df_list), 1)){
  df_spots <- rbind(df_spots, data.frame(df_list[[i]], id = i))
}

# Convert Spots df to numeric
df_spots <- df_spots[, c("id","TRACK_ID", "POSITION_X", "POSITION_Y", "POSITION_T")]
df_spots$POSITION_X <- as.numeric(df_spots$POSITION_X)
df_spots$POSITION_Y <- as.numeric(df_spots$POSITION_Y)
df_spots$POSITION_T <- as.numeric(df_spots$POSITION_T)
df_spots <- na.omit(df_spots)

# Assign n_id
df_tracks$n_id <- paste(df_tracks$id, df_tracks$TRACK_INDEX)
df_tracks$case <- paste(df_tracks$Genotype, df_tracks$Stim)
df_spots$n_id <- paste(df_spots$id, df_spots$TRACK_ID)

# Merge spot and track dataset to combine average speed and position in one dataframe
spot_track <- merge(df_spots, df_tracks, by = "n_id")
spot_track <- spot_track[, c("n_id", "POSITION_X", "POSITION_Y", "POSITION_T", "TRACK_MEAN_SPEED", "case", "Stim")]

# Pixel in um und pixel/10min in um/h umrechnen
spot_track$TRACK_MEAN_SPEED <- spot_track$TRACK_MEAN_SPEED*14.781605114
spot_track$POSITION_X <- spot_track$POSITION_X*2.46360085
spot_track$POSITION_Y <- spot_track$POSITION_Y*2.46360085

# Extract starting point from 
start_point <- spot_track %>% group_by(n_id) %>% dplyr::slice(which.min(spot_track$POSITION_T))
names(start_point)[names(start_point) == "POSITION_X"] <- "start_x"
names(start_point)[names(start_point) == "POSITION_Y"] <- "start_y"
start_point <- start_point[,c("n_id", "start_x", "start_y")]

# Bring Start point back to spot_track
spot_track <- merge(spot_track, start_point, by = c("n_id"))
spot_track <- spot_track[order(spot_track$n_id, spot_track$POSITION_T),]

# Subtract starting point from each coordinate
spot_track$zero_x <- spot_track$POSITION_X-spot_track$start_x
spot_track$zero_y <- spot_track$POSITION_Y-spot_track$start_y

# Rename Groups
spot_track <- spot_track[spot_track$Stim == "CT",]

spot_track$group <- 
  mgsub(spot_track$case,
        c("B2 CT", "B3 CT", "B4 CT"), 
        c("KO", "WT", "Inact."))

# Order by timepoint and group
spot_track <- transform(spot_track,
                        group=factor(group,levels=c("KO", 
                                                    "WT", 
                                                    "Inact.")))

########################## Fig. 6C - Trajectory Maps ###########################

ggplot(spot_track, aes(x = spot_track$zero_x, y = spot_track$zero_y, 
                       group = spot_track$n_id, color = spot_track$TRACK_MEAN_SPEED))+
  geom_hline(aes(yintercept = 0), size = 0.3)+
  geom_vline(aes(xintercept = 0), size = 0.3)+
  geom_path(size = 0.3)+    
  xlim(-350, 350)+
  ylim(-350, 350)+
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
        axis.title.y = element_text(size =8, face= "plain", family = "Arial", angle = 90, margin = margin(r = 1)),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", angle = 0),
        legend.text = element_text(size =6, face= "plain", family = "Arial"),
        legend.title = element_text(size =8, face= "plain", family = "Arial"),
        legend.key.height = unit(0.1, "inches"), #change legend key height
        legend.key.width = unit(0.05, "inches"),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_line(colour = "grey", size = 0.05))+
  facet_grid(. ~group)+
  theme(strip.text.y = element_blank())
ggsave(paste0("Output/Fig. 4J - Trajectory Maps.tiff"), width = 3.2, height = 1.1, units = "in", dpi = 1000)

############################ Fig. 6C - Analysis ################################

# Assign position well
df_tracks$well <- paste(df_tracks$case, df_tracks$id)

# Weight speed of each track by th captured spots in the track
df_tracks$speedxedges <- df_tracks$NUMBER_EDGES*df_tracks$TRACK_MEAN_SPEED
Erg <- aggregate(df_tracks[, c("NUMBER_EDGES", "speedxedges")], list(df_tracks$well), mean)
Erg$average_speed_per_frame <- Erg$speedxedges/Erg$NUMBER_EDGES
writexl::write_xlsx(path = "Output/Fig. 4J - Analysis.xlsx", x = Erg)
