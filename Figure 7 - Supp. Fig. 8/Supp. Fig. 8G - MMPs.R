library(ggplot2)
setwd("/Users/larskoch/Library/CloudStorage/Dropbox/3. Adamts12/Revision/10) MMP,ADAM,ADAMTS in massspec_bulk_visium")

# Read Dataset
pro <- readxl::read_excel("230531 log2 transformed intensities proteome two sample t test.xlsx")

# Plot all ADAMTS & ADAM & MMP
MMP <- pro[grep("^(ADAM|MMP|ADAMTS)", pro$PG.Genes),]

# Supp. Fig. 8G
ggplot(MMP, aes(x = -MMP$`Student's T-test Difference NT_ADAMTS12ko`,
                y = reorder(MMP$PG.Genes, MMP$`Student's T-test Difference NT_ADAMTS12ko`)))+ 
  geom_bar(stat = "identity", aes(fill =-MMP$`Student's T-test Difference NT_ADAMTS12ko`))+
  geom_vline(xintercept = 0,
             color = "black", size = 0.3)+
  labs(title = "DEP") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  scale_x_continuous(breaks = c(-4,-2,0,2))+
  theme_minimal() +
  theme(plot.title = element_text(size = 9, face="plain", hjust = 0.5,),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 8, face = "plain", family="Arial"),
        axis.text.x = element_text(size =6, hjust = 1, face= "plain", family= "Arial", color = "black"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial", color = "black"),
        axis.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_blank() ,
        legend.key.width= unit(0.4, 'cm'),)+
  xlab("t-value")
ggsave("Output/ADAMTS_ADAM_MMP Expression in Mass Spect.tiff", width = 1.4, height = 2.4, units = "in", bg = "white")
