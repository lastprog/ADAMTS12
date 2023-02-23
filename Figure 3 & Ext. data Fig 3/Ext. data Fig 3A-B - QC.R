# Written in 2023 by Lars Koch
# Contact: https://www.linkedin.com/in/lars-koch-93844023a/
# Publication: "ADAMTS12 promotes fibrosis by restructuring ECM to enable activation 
# and migration of injury-responsive fibroblasts"

library(tidyverse)
library(ggplot2)
library(mgsub)
library(stringr)
library(ggpubr)
setwd("xxx") # Set your working directory here

# Load metrics summary data
samples <- c("KH42", "KH43", "KH44", "KH45")
metric_sum <- list()
for(i in samples){
    metric_sum[[i]] <- read.csv(paste0("Input/metrics_summary ", i, ".csv"))
}

# Create comprehensive dataframe of metrics data of all samples
df <- data.frame(KH42 = t(metric_sum[["KH42"]]),
                 KH43 = t(metric_sum[["KH43"]]),
                 KH44 = t(metric_sum[["KH44"]]),
                 KH45 = t(metric_sum[["KH45"]]))

# Prepare plotting
df_long <- gather(df, key = "sample", value = "value")
df_long$rownames <- rep(rownames(df), length(samples))
# Rename KH42, KH43 etc to "Adamts12-/- #1" and "Adamts12-/- #2"
df_long$sample <- mgsub(df_long$sample, c("KH42", "KH43"), c("Adamts12-/- #1", "Adamts12-/- #2"))
df_long$sample <- mgsub(df_long$sample, c("KH44", "KH45"), c("WT #1", "WT #2"))
df_long$sample <- factor(df_long$sample, levels = c("Adamts12-/- #2", "Adamts12-/- #1", "WT #2", "WT #1"))
df_long$group <- str_split_fixed(df_long$sample, " ", 2)[,1]
df_long$value <- as.numeric(df_long$value)

################### Ext. Data Fig 3B - Total Detected Genes ###################

to_ge_det <- df_long[df_long$rownames == "Total.Genes.Detected",]
p1 <- ggplot(to_ge_det, aes(x = to_ge_det$sample, y = to_ge_det$value/1000, fill = to_ge_det$group))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = c(0,0), breaks = c(0,5,10,15))+
  ggtitle("Total Detected Genes")+
  ylab("[in 1000]")+
  theme_classic()+
  theme(plot.title = element_text(size =9, angle = 0, hjust = 0.5, face= "plain", family= "Arial", color = "black"), 
        axis.title = element_blank(),
        axis.title.x = element_text(size =8,  face= "plain", family= "Arial", color = "black"),
        axis.text.x = element_text(size =6, face= "plain", family= "Arial", color = "black"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial", color = "black"),
        legend.position = "none")+
  coord_flip()

######################## Ext. Data Fig 3C - Valid UMIs ########################

val_umis <- df_long[df_long$rownames == "Valid.UMIs",]
p2 <- ggplot(val_umis, aes(x = val_umis$sample, y = val_umis$value*100, fill = val_umis$group))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,50,100))+
  ggtitle("Valid UMIs")+
  ylab("Valid UMIs [%]")+
  theme_classic()+
  theme(plot.title = element_text(size =9, angle = 0, hjust = 0.5, face= "plain", family= "Arial", color = "black"), 
        axis.title = element_blank(),
        axis.title.x = element_text(size =8,  face= "plain", family= "Arial", color = "black"),
        axis.text.x = element_text(size =6,  face= "plain", family= "Arial", color = "black"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial", color = "black"),
        legend.position = "none")+
  coord_flip()

# Arrange and save plots
ggarrange(p1, p2, ncol = 1)
ggsave(filename = "Output/Supp. Fig 3C QC.svg",height = 2.2, width = 1.6, units = "in")







