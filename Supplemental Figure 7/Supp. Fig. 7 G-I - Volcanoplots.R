# Written in 2023 by Lars Koch
# Contact: https://www.linkedin.com/in/lars-koch-93844023a/
# Publication: "ADAMTS12 promotes fibrosis by restructuring ECM to enable activation 
# and migration of injury-responsive fibroblasts"

library(ggplot2)
library(ggrepel)
setwd("xxx") # Set your working directory here

# Vulcano_nice function with settable y-axis limit, because -log(P-Val) = Infinity
volcano_nice <- function (df, hAss = 0.05, FCIndex, pValIndex, IDIndex, vAss = NULL,
                                      label = FALSE, straight = FALSE, nlabels, manual_labels = NA, ylimab = 20, xlimab = 5)
{
  df <- df[complete.cases(df), ]
  names(df)[1] <- "X1"
  hAssOri <- hAss
  hAss <- -log(hAss)
  names(df) <- gsub("P.Val", "FDR", names(df))
  names(df)[FCIndex] <- "logFC"
  names(df)[pValIndex] <- "P.Val"
  if (max(abs(df[, FCIndex])) >= 1) {
    xlimAbs <- xlimab
    ylimAbs <- ylimab
  }
  else {
    xlimAbs <- xlimab
    ylimAbs <- ylimab
  }
  if (is.null(vAss)) {
    vAss <- xlimAbs/10
  }
  xneg <- function(x) abs(hAss - 1 + x/(x + vAss))
  xpos <- function(x) abs(hAss - 1 + x/(x - vAss))
  test <- function(x, y, vAss) {
    if (x < -vAss) {
      if (xneg(x) < -log(y)) {
        return("1")
      }
      else {
        return("0")
      }
    }
    else {
      if (x > vAss) {
        if (xpos(x) < -log(y)) {
          return("1")
        }
        else {
          return("0")
        }
      }
      else {
        return("0")
      }
    }
  }
  if (straight) {
    df$couleur <- ifelse(abs(df$logFC) >= vAss & df$P.Val <=
                           hAssOri, "1", "0")
  }
  else {
    df$couleur <- "0"
    df$couleur <- apply(df, 1, FUN = function(x) test(as.numeric(x[FCIndex]),
                                                      as.numeric(x[pValIndex]), vAss))
  }
  df <- df[order(df$P.Val, decreasing = F), ]
  df$condLabel <- df[, IDIndex]
  df[df$couleur == "0", "condLabel"] <- NA
  labels_to_keep <- c(df[c(1:nlabels), "condLabel"],manual_labels)
  df[!(df$condLabel %in% labels_to_keep), "condLabel"] <- NA
  df$couleur <- ifelse(df$couleur == "1" & df$logFC < 0, "2",
                       df$couleur)
  if (label) {
    a <- ggplot(df, aes(x = logFC, y = -log(P.Val), color = couleur)) +
      geom_point(alpha = 1, size = 0.1) + geom_label_repel(aes(label = condLabel)) +
      stat_function(fun = xneg, xlim = c(-xlimAbs, -vAss),
                    color = "black", alpha = 0.7) + ylim(c(0, ylimAbs)) +
      xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
      scale_colour_manual(values = c("grey30", "indianred",
                                     "darkblue")) + theme_minimal() + theme(legend.position = "none")+
      theme(axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
            axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"))
  }
  else {
    if (straight) {
      a <- ggplot(df, aes(x = logFC, y = -log(P.Val),
                          color = couleur)) + geom_point(alpha = 1, size = 0.1) +
        geom_vline(xintercept = -vAss, color = "blue") +
        geom_vline(xintercept = vAss, color = "blue") +
        ylim(c(0, ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) +
        geom_hline(yintercept = hAss, color = "red") +
        scale_colour_manual(values = c("grey30", "indianred",
                                       "darkblue")) + theme_minimal() + theme(legend.position = "none")+
        theme(axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
              axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"))
    }
    else {
      a <- ggplot(df, aes(x = logFC, y = -log(P.Val),
                          color = couleur)) + geom_point(alpha = 1, size = 0.1) +
        stat_function(fun = xneg, xlim = c(-xlimAbs,
                                           -vAss), color = "black", alpha = 0.7) + ylim(c(0,
                                                                                          ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                                                                                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
        scale_colour_manual(values = c("grey30", "indianred",
                                       "darkblue")) + theme_minimal() + theme(legend.position = "none")+
        theme(axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
              axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"))
    }
  }
  return(a)
}

######################### Ext data Fig 4l - Act vs KO ##########################

txt <- "WTvsKO"

# Read Data
deg <- read.table(paste0("Input/ADAMTS12_OvExp_00_DESEQ2_", txt,".txt"), header = T)
colnames(deg) <- c("ID", "base.Meanx", "log2FC", "lfc.SE", "t", "pval", "adj.pval", "index")
write.csv(deg, paste0("Supp. table 5 - DEGs Bulk Seq", txt,".csv"))

# For annotation
volcano_nice(df = deg, hAss = 0.05, FCIndex = 3,pValIndex = 6, vAss = 0.3,
             IDIndex = 1, label = T,nlabels = 2000, ylimab = 30, xlimab = 5)+
  ggtitle(txt)+
  labs(x="logFC", y="-log(P-Val.)")+
  theme(plot.title = element_text(size = 9, family ="Arial", hjust = 0.5),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_line(colour = "grey", size = 0.05))
ggsave(filename = paste0("Output/Ext data Fig 4 l_",txt,".jpg"), width = 30, height =30)

# For plotting
volcano_nice(df = deg, hAss = 0.05, FCIndex = 3,pValIndex = 6, vAss = 0.3,
             IDIndex = 1, label = F,nlabels = 20, ylimab = 30, xlimab = 5)+
  ggtitle(txt)+
  labs(x="logFC", y="-log(P-Val.)")+
  theme(plot.title = element_text(size = 9, family ="Arial", hjust = 0.5),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_line(colour = "grey", size = 0.05))
ggsave(filename = paste0("Output/Ext data Fig 4 l_",txt,".svg"), width = 2, height = 2)

######################## Ext data Fig 4l - Act vs Inact ########################

txt <- "WTvsInact"

# Read Data
deg <- read.table(paste0("Input/ADAMTS12_OvExp_00_DESEQ2_", txt,".txt"), header = T)
colnames(deg)<-c("ID", "base.Meanx", "log2FC", "lfc.SE", "t", "pval", "adj.pval", "index")
write.csv(deg, paste0("Supp. table 5 - DEGs Bulk Seq", txt,".csv"))

# For annotation
volcano_nice(df = deg, hAss = 0.05, FCIndex = 3,pValIndex = 6, vAss = 0.3,
             IDIndex = 1, label = T,nlabels = 1000, ylimab = 31, xlimab = 5)+
  ggtitle(txt)+
  labs(x="logFC", y="-log(P-Val.)")+
  theme(plot.title = element_text(size = 9, family ="Arial", hjust = 0.5),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_line(colour = "grey", size = 0.05))
ggsave(filename = paste0("Output/Ext data Fig 4 l_",txt,".jpg"), width = 50, height =50, limitsize = F)

# For plotting
volcano_nice(df = deg, hAss = 0.05, FCIndex = 3,pValIndex = 6, vAss = 0.3,
             IDIndex = 1, label = F,nlabels = 200, ylimab = 31, xlimab = 5)+
  ggtitle(txt)+
  labs(x="logFC", y="-log(P-Val.)")+
  theme(plot.title = element_text(size = 9, family ="Arial", hjust = 0.5),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_line(colour = "grey", size = 0.05))
ggsave(filename = paste0("Output/Ext data Fig 4 l_",txt,".svg"), width = 2, height = 2)

######################## Ext data Fig 4l - Act vs Inact ########################

txt <- "InactvsKO"

# Read Data
deg <- read.table(paste0("Input/ADAMTS12_OvExp_00_DESEQ2_", txt,".txt"), header = T)
colnames(deg)<-c("ID", "base.Meanx", "log2FC", "lfc.SE", "t", "pval", "adj.pval", "index")
write.csv(deg, paste0("Supp. table 5 - DEGs Bulk Seq", txt,".csv"))

# For annotation
volcano_nice(df = deg, hAss = 0.05, FCIndex = 3,pValIndex = 6, vAss = 0.3,
           IDIndex = 1, label = T,nlabels = 200, ylimab = 25, xlimab = 5, manual_labels = c("FOS", "JUNB", "ATF3"))+
  ggtitle(txt)+
  labs(x="logFC", y="-log(P-Val.)")+
  theme(plot.title = element_text(size = 9, family ="Arial", hjust = 0.5),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_line(colour = "grey", size = 0.05))
ggsave(filename = paste0("Output/Ext data Fig 4 l_",txt,".jpg"), width = 30, height =30)
  
# For plotting
volcano_nice(df = deg, hAss = 0.05, FCIndex = 3,pValIndex = 6, vAss = 0.3,
             IDIndex = 1, label = F,nlabels = 200, ylimab = 25, xlimab = 5)+
  ggtitle(txt)+
  labs(x="logFC", y="-log(P-Val.)")+
  theme(plot.title = element_text(size = 9, family ="Arial", hjust = 0.5),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_line(colour = "grey", size = 0.05))
ggsave(filename = paste0("Output/Ext data Fig 4 l_",txt,".svg"), width = 2, height = 2)

