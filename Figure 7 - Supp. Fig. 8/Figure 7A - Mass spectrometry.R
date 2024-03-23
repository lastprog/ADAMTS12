# Written in 2023 by Lars Koch
# Contact: https://www.linkedin.com/in/lars-koch-93844023a/
# Publication: "ADAMTS12 promotes fibrosis by restructuring ECM to enable activation 
# and migration of injury-responsive fibroblasts"

library(readxl)
library(ggplot2)
setwd("xxx") # Set your working directory here

df <- readxl::read_excel("Input/220619_DEG_ProteinLevel.xlsx", col_names = T)
colnames(df) <- df[1,]
df <- df[-1,]

# Arrange Dataframe
df$`Student's T-test Difference NT_ADAMTS12ko` <- as.numeric(gsub(",", ".", df$`Student's T-test Difference NT_ADAMTS12ko`))
df$`Student's T-test Difference NT_ADAMTS12ko` <- df$`Student's T-test Difference NT_ADAMTS12ko`*-1
df$`Student's T-test q-value NT_ADAMTS12ko` <- as.numeric(gsub(",", ".", df$`Student's T-test q-value NT_ADAMTS12ko`))
colnames(df) <- c("significant", "q-value", "Diff_KO_vs_WT", "Protein_groups", "Genes", "UniProtIds")
df <- df[order(df$Diff_KO_vs_WT),]

write.csv(df, "Output/Supp. table 4 - DEP mass spectrometry.csv")

## Top 10 up and down
top10 <- rbind(head(df, 10), tail(df,10))

# Figure 7A
ggplot(top10, aes(x = reorder(top10$Genes, -top10$Diff_KO_vs_WT),
                  y = top10$Diff_KO_vs_WT)) + 
  geom_bar(stat = "identity", aes(fill =top10$Diff_KO_vs_WT))+
  geom_hline(yintercept = 0,
             color = "black", size = 0.3)+
  labs(title = "DEP") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  scale_y_continuous(breaks = c(-6,-3,0,3))+
  theme_minimal() +
  theme(plot.title = element_text(size = 9, face="plain", hjust = 0.5,),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 8, face = "plain", family="Arial"),
        axis.text.x = element_text(size =6,, hjust = 1, face= "plain", family= "Arial", color = "black"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial", color = "black"),
        axis.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_blank() ,
        legend.key.width= unit(0.4, 'cm'),) +
  ylab("Diff. expression")+
  coord_flip()
ggsave(filename = "C:/Users/Lars/Dropbox/3. Adamts12/Figures/Figures (ai)/Figure 4/Figure 4 Daten & Abbildungen/mass_spectrometry_DEP.svg",  
       width = 1.4, height = 2.4, units = "in")

