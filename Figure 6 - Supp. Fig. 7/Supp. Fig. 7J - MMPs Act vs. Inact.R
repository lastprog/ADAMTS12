library(ggplot2)

# Read dataset
deg_act_inact <- read.csv(
  "/Users/larskoch/Library/CloudStorage/Dropbox/3. Adamts12/Figures/Figures (ai)/Supplemental Tables/Supp. table 6 - DEGs Bulk Seq WT_Veh.vs.Inact_Veh.csv")

# Select all ADAMTS & ADAM & MMP
MMP <- deg_act_inact[grep("^(ADAM|MMP|ADAMTS)\\d+", deg_act_inact$ID),]

# Supp. Fig. 7J
ggplot(MMP, aes(x = -MMP$t,
                y = reorder(MMP$ID, -MMP$t)))+ 
  geom_bar(stat = "identity", aes(fill =-MMP$t))+
  geom_vline(xintercept = 0,
             color = "black", size = 0.3)+
  labs(title = "Metalloproteases") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal()+
  theme(plot.title = element_text(size = 9, face="plain", hjust = 0.5,),
        axis.title = element_blank(),
        axis.title.y = element_text(angle = 90, size = 8, face = "plain", family="Arial"),
        axis.text.x = element_text(angle = 45, size =6, hjust = 1, face= "plain", family= "Arial", color = "black"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial", color = "black"),
        axis.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_blank() ,
        legend.key.width= unit(0.4, 'cm'),)+
  xlab("t-value")+
  coord_flip()
ggsave("/Users/larskoch/Library/CloudStorage/Dropbox/3. Adamts12/Revision/0) Figures, Tables and Western Blots/Figure 6 Daten und Abbildungen/ADAMTS_ADAM_MMP Expression in bulk act_inact.svg", width = 4.5, height = 1.4, units = "in", bg = "white")
