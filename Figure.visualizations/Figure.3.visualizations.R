library(dplyr)
library(ggplot2)
library(stringr)
library(scales)
library(CellChat)
library(tidyverse)

#Read in MASC results of DRE vs DCE
setwd("~")
MASCinput <- read_excel("Supplementary table 2.xlsx", sheet = 3)
MASCinput$cluster[MASCinput$cluster == "NK"] <- "NKCD56dim"
MASCinput <- MASCinput %>% filter(size > 100, CI.upper < 16, cluster != "HSPC")
MASCinput$cluster <- factor(MASCinput$cluster,
                            levels = c("B", "Bnaive", "Bintermediate", "Bmemory",
                                       "CD4T", "CD4Naive", "CD4TCM", "CD4TEM", "CD4CTL",
                                       "Treg", "CD8T", "CD8Naive", "CD8TCM", "CD8TEM",
                                       "otherT", "gdT", "MAIT", "Mono", "CD14Mono",
                                       "CD16Mono", "DC", "cDC2", "pDC", "NK", "NKCD56dim"))
MASCinput$Significance <- "Not significant"
MASCinput$Significance[MASCinput$CI.lower > 1] <- "Enriched"
MASCinput$Significance[MASCinput$CI.upper < 1] <- "Reduced"
MASCinput$Significance <- factor(MASCinput$Significance,
                              levels = c("Enriched", "Not significant", "Reduced"))

ggplot(MASCinput %>% filter(!cluster %in% c("B", "CD4T", "CD8T", "otherT", "Mono", "DC", "NK")),
       aes(x = cluster, y = log2(OR), ymin = log2(CI.lower),
                               ymax = log2(CI.upper), col = Significance)) +
  geom_pointrange() +
  xlab("Cell type") + ylab("log2(Odds Ratio)") + ylim(c(-1.25,1)) +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggtitle("Enrichment of immune cell types in DRE compared to DCE") +
  scale_color_manual(values = c("red", "black", "blue"))

#GSEA plot of DRE vs DCE
sheet_names <- excel_sheets("Supplementary table 9.xlsx")
gsealist <- map(sheet_names,
                       ~ read_excel("Supplementary table 9.xlsx", sheet = .x)) %>% 
  set_names(sheet_names)
L1gseaoutput <- gsealist[c(1:7)]
L2gseaoutput <- gsealist[c(8:25)]

fullgsea <- bind_rows(L1gseaoutput$CD8T, L1gseaoutput$DC, L1gseaoutput$Mono,
                      L2gseaoutput$Bintermediate, L2gseaoutput$CD14Mono,
                      L2gseaoutput$CD16Mono, L2gseaoutput$CD4Naive,
                      L2gseaoutput$CD4CTL, L2gseaoutput$CD8Naive,
                      L2gseaoutput$pDC)

pathways <- fullgsea$Description[c(5,11,12,20,21,27,28,32,39,48,50,72,111,126,146,328,335)]

ggplot(fullgsea %>% filter(Description %in% pathways),
       aes(x = CellType, y = Description, size = p.adjust, col = NES)) +
  geom_point() + theme_bw() + ylab(label = " ") +
  ggtitle("Selected pathways in DRE vs DCE") +
  scale_y_discrete(labels = wrap_format(40))+
  theme(plot.title = element_text(hjust = 0.75)) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(axis.title.x = element_blank()) +
  scale_color_gradient(low = "#132B43", high = "#3670A0")

#CellChat of PCE vs WCE
load("Cellchat.new.Seurat.object.l2.DCE.RData")
load("Cellchat.new.Seurat.object.l2.DRE.RData")
load("cellchat.new.l2.DRE.DCE.RData")
cellchatl2cr <- liftCellChat(cellchatl2cr)

signaling <- c("EPHB", "COMPLEMENT", "IGF", "MHC-I", "CCL", "ICAM",
               "GALECTIN", "SELL", "SELPLG", "IFN-II", "CD80", "CD86", "IL1", "MHC-II",
               "IL6", "IL10")
ranknetwcpc <- rankNet(cellchatl2cr, mode = "comparison", measure = "weight", do.stat = TRUE,
                       signaling = signaling,
                       color.use = c("#FFA500", "#7EC0EE"), stacked = T,
                       title = "Selected signaling in DCE vs DRE")
ranknetwcpc$theme$axis.text.y$colour <- c(rep("#7EC0EE", 13), rep("#FFA500", 3))

sources <- levels(cellchat2@idents)[c(2:8,10:12,14:16,18,21,22,24,25,28)]
netVisual_aggregate(cellchat1, signaling = "IL10", layout = "circle",
                    remove.isolate = T)
netVisual_aggregate(cellchat2, signaling = "IL10", layout = "circle",
                    remove.isolate = T)
