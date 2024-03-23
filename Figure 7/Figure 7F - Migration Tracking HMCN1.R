library(stringr)
library(dplyr)
library(ggplot2)
library(mgsub)
library(xml2)
library(tibble)
setwd("/Users/larskoch/Library/CloudStorage/Dropbox/3. Adamts12/Revision/HNCM1 siRNA migration tracking/Migration Tracking HMCN1 siRNA V4")

# Import and arrange track datafrmae
dir <- "Input/csv"
files <- list.files(dir)
df_list <- lapply(paste0(dir, "/", files), function(x) read.csv(x))

# Convert list into dataframe, with filename as ID
df_tracks <- data.frame()
for(i in seq(1, length(df_list), 1)){
  df_tracks <- rbind(df_tracks, data.frame(df_list[[i]], id = i))
}

files <- sub(".csv", " ", files)
files <- as.numeric(files)
df <- data.frame(id = 1:length(files),
                 files = files)
df_tracks <- merge(df_tracks, df, "id")
df_tracks <- df_tracks[, c("id","files", "TRACK_ID", "TRACK_MEAN_SPEED", "NUMBER_SPOTS")]

# Assign group with multipoints
doc <- read_xml("Input/xml names/Migration Tracking all.xml")

names <- data.frame()
for(i in 0:355){
  point_name <- paste0(".//", sprintf("Point%05d",i))
  temp <- xml_find_all(doc, point_name)
  name <- xml_attrs(xml_child(temp[[1]], 2))[["value"]]
  names <- rbind(names, c(name, i+1))
}
colnames(names) <- c("name", "id")

df_tracks <- merge(df_tracks, names, by = "id")

# Rename
df_tracks$genotype <- str_split_fixed(df_tracks$name, " ", 4)[,1]
df_tracks$stim <- str_split_fixed(df_tracks$name, " ", 4)[,2]
df_tracks$well <- str_split_fixed(df_tracks$name, " ", 4)[,3]
df_tracks$group <- paste(df_tracks$genotype, df_tracks$stim, df_tracks$well)
df_tracks$n_id <- paste(df_tracks$id, df_tracks$TRACK_ID)

# Mean
df_tracks$NUMBER_EDGES <- df_tracks$NUMBER_SPOTS -1
df_tracks$speedxedges <- df_tracks$NUMBER_EDGES*df_tracks$TRACK_MEAN_SPEED

Erg <- aggregate(df_tracks[, c("NUMBER_EDGES", "speedxedges")], list(df_tracks$group, df_tracks$files), mean)
Erg$average_speed_per_frame <- Erg$speedxedges/Erg$NUMBER_EDGES

# Calculate speed in um/h (from um/s)
Erg$speed_um_per_h <- Erg$average_speed_per_frame*661.5997 # per 10 minutes (frame)
Erg$speed_um_per_h <- Erg$speed_um_per_h*6 # per hour

# Rename
Erg$genotype <- str_split_fixed(Erg$Group.1, " ",3)[,1]
Erg$condition <- str_split_fixed(Erg$Group.1, " ",3)[,2]
writexl::write_xlsx(path = "Output/B2-B3 HNCM1 V5 speed_um_per_h_publish.xlsx", x = Erg)

# Anova
Erg$group <- paste(Erg$genotype, Erg$condition)
anova <- TukeyHSD(aov(speed_um_per_h ~ group, data = Erg))
av <- as.data.frame(anova[[1]])
av <- rownames_to_column(av)
writexl::write_xlsx(path = "Output/B2-B3 HNCM1 V5 reanalysis num_spots_av.xlsx", av)

##################### Start Calculating Trajectory maps #######################
df_spots <- read.csv("/Users/larskoch/Library/CloudStorage/Dropbox/3. Adamts12/Revision/HNCM1 siRNA migration tracking/Migration Tracking HMCN1 siRNA V4/Output/df_spots all V4.csv")

spot_track <- merge(df_spots, df_tracks, by = "n_id")
spot_track$genotype <- str_split_fixed(spot_track$name, " ", 4)[,1]
spot_track$stim <- str_split_fixed(spot_track$name, " ", 4)[,2]

spot_track <- spot_track[, c("n_id", "POSITION_X", "POSITION_Y", "POSITION_T", "TRACK_MEAN_SPEED", "genotype", "stim")]

# Pixel/10min to um/h
spot_track$TRACK_MEAN_SPEED <- spot_track$TRACK_MEAN_SPEED*661.5997 # Per 10 minutes (frame)
spot_track$TRACK_MEAN_SPEED <- spot_track$TRACK_MEAN_SPEED*6 # Per hour

# Extract starting point from each track
spot_track <- group_by(spot_track, spot_track$n_id)

spot_track$POSITION_T <- as.numeric(as.character(spot_track$POSITION_T))
spot_track$POSITION_X <- as.numeric(as.character(spot_track$POSITION_X))
spot_track$POSITION_Y <- as.numeric(as.character(spot_track$POSITION_Y))
spot_track <- spot_track[,c("n_id", "POSITION_X", "POSITION_Y", "POSITION_T", "TRACK_MEAN_SPEED", "genotype", "stim")]

spot_track$fov <- as.numeric(as.character(str_split_fixed(spot_track$n_id, " ", 2)[,1]))
spot_track$track_id <- as.numeric(as.character(str_split_fixed(spot_track$n_id, " ", 2)[,2]))
start_point <- data.frame()
for(i in unique(spot_track$fov)){
  tracks <- pull(unique(spot_track[spot_track$fov == i,"track_id"]))
  for(t in tracks){
    df_subs <- spot_track[spot_track$fov == i & spot_track$track_id == t, ]
    min <- which.min(df_subs$POSITION_T)
    start_point <- rbind(start_point, df_subs[min,])
  }
  print(paste(i, "done", Sys.time()))
}

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
spot_track$group <- paste(spot_track$genotype, spot_track$stim)

# Order by timepoint and group
spot_track <- transform(spot_track,
                        group=factor(group,levels=c("B2 CT", 
                                                    "B2 siRNA", 
                                                    "B3 CT",
                                                    "B3 siRNA")))

# Figure 7F
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
        axis.title.x = element_text(size =8, face= "plain", family = "Arial"),
        axis.title.y = element_text(size =8, face= "plain", family = "Arial", angle = 90),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", angle = 0),
        legend.text = element_text(size =6, face= "plain", family = "Arial"),
        legend.title = element_text(size =8, face= "plain", family = "Arial"),
        legend.key.height = unit(0.1, "inches"), 
        legend.key.width = unit(0.05, "inches"),
        panel.grid.major = element_line(colour = "grey", size = 0.1),
        panel.grid.minor = element_line(colour = "grey", size = 0.1))+
  facet_grid(genotype ~ stim)

ggsave("Output/trajectory maps_siRNA_03.tiff", width = 2.7, height = 2.2, units = "in", dpi = 1000)




