# Written in 2023 by Lars Koch
# Contact: https://www.linkedin.com/in/lars-koch-93844023a/
# Publication: "ADAMTS12 promotes fibrosis by restructuring ECM to enable activation 
# and migration of injury-responsive fibroblasts"

library(readxl)
library(matrixStats)
library(stringr)
library(ggplot2)
setwd("xxx") # Set your working directory here

# Read and arrange PDGFRb probabilities predicted using ilastik
pdgfrb <- read_xlsx("Input/PDGFRb Probs.xlsx", col_names = F)
coln <- paste(pdgfrb[1,], pdgfrb[2,], pdgfrb[3,])
pdgfrb <- pdgfrb[-c(1:4),]
colnames(pdgfrb) <- coln

# Read and arrange ADAMTS12 probabilities predicted using ilastik
adamts12 <- read_xlsx("Input/ADAMTS12 Probs.xlsx", col_names = F)
adamts12 <- adamts12[-c(1:4),]
colnames(adamts12) <- coln

# Create Boolean from probability, probability >50% means positive cell
pdgfrb_b <- ifelse(pdgfrb > 0.5, T, F)
adamts12_b <- ifelse(adamts12 > 0.5, T, F)

# Shared positivity
pdgfrb_adamts12_pos <- as.matrix(ifelse(pdgfrb_b == T & adamts12_b == T, T, F))
pdgfrb_neg_adamts12_pos <- as.matrix(ifelse(pdgfrb_b == F & adamts12_b == T, T, F))

# Count occurence 
counts <- data.frame(row = str_split_fixed(coln, " ", 3)[,1],
                     col = str_split_fixed(coln, " ", 3)[,2],
                     pic = str_split_fixed(coln, " ", 3)[,3],
                     pp = colCounts(pdgfrb_b, value = T, na.rm = T),
                     pn = colCounts(adamts12_b, value = F, na.rm = T),
                     pp_ap = colCounts(pdgfrb_adamts12_pos, value = T, na.rm = T),
                     pn_ap = colCounts(pdgfrb_neg_adamts12_pos, value = T, na.rm = T))

# Shared positivity in percent of PDGFRb+ and PDGFRb-
counts$pp_ap_dev_pp <- counts$pp_ap/counts$pp
counts$pn_ap_dev_pn <- counts$pp_ap/counts$pn

# Summarize per kidney
mean_per_kid <- aggregate(pp_ap_dev_pp ~ row + col, counts, FUN = mean)
mean_per_kid$pp_ap_dev_pn <- aggregate(pn_ap_dev_pn ~ row + col, counts, FUN = mean)[,"pn_ap_dev_pn"]
mean_per_kid$pp_ap_dev_pp <- mean_per_kid$pp_ap_dev_pp*100
mean_per_kid$pp_ap_dev_pn <- mean_per_kid$pp_ap_dev_pn*100

# Exclude kidneys with tumor
mean_per_kid$kid <- paste(mean_per_kid$row, mean_per_kid$col)
ex <- c("2 4", "2 6", "4 1", "4 4", "6 1", "6 2", "6 6")
mean_per_kid <- mean_per_kid[!(mean_per_kid$kid %in% ex),]

writexl::write_xlsx(mean_per_kid, path = "Output/Occurence by pn pp per kid.xlsx")

# Calculate T-Test
t.test(x = mean_per_kid$pp_ap_dev_pp, y = mean_per_kid$pp_ap_dev_pn, paired = T)


################# Fig 1G - ADAMTS12 occurence in PDGFRb+/- cells ###############

mean_longer <- data.frame(name = c(rep("PDGFRb+", nrow(mean_per_kid)),rep("PDGFRb-", nrow(mean_per_kid))),
                          percent = c(mean_per_kid$pp_ap_dev_pp, mean_per_kid$pp_ap_dev_pn),
                          id= c(seq(1, nrow(mean_per_kid), 1), seq(1, nrow(mean_per_kid), 1)))
mean_longer$name <- factor(mean_longer$name, levels = c("PDGFRb+", "PDGFRb-"))

# Plot
ggplot(mean_longer)+
  geom_boxplot(aes(y = mean_longer$percent, x = mean_longer$name), outlier.size =  0.3)+
  geom_line(aes(y = mean_longer$percent, x = mean_longer$name, group = mean_longer$id, color = mean_longer$name),color = "black", lwd = 0.01)+
  geom_point(aes(y = mean_longer$percent, x = mean_longer$name, group = mean_longer$id, color = mean_longer$name),size = 0.3)+
  ylab("ADAMTS12 expr. cells [%]")+
  scale_y_continuous(expand = c(0, 0),)+
  scale_color_manual(values = c("PDGFRb-" = "darkblue", "PDGFRb+" = "indianred"))+
  theme_classic()+
  theme(plot.title = element_blank(), 
        axis.title.x = element_text(size = 8, face = "plain", family="Arial"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size =6, angle = 0, hjust = 1, face= "plain", family= "Arial", color = "black"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial", color = "black"),
        axis.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.position = "none",
        legend.key.width= unit(0.4, 'cm'))+
  coord_flip()
ggsave("Output/PDGFRb pos neg boxplot_dots and lines.svg", height = 1, width = 2)  




