library(ggplot2)
library(shiny)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(Seurat)
library(ggrepel)
library(DT)
library(ggpubr)
library(dplyr)
library(gridExtra)
library(visNetwork)
library(bslib)
library(stringr)
library(CellChat)
library(cowplot)
library(scRepertoire)
library(ggiraph)
library(reticulate)

#############

load(file = "shiny_objects_3_1231.Rdata")
load(file="MASC_plots_tables.RData")
markergenes <- rownames(data_matrix)
load("New.TCR.cloneType.data.from.Seurat.object.RData")
load("Cellchat.new.Seurat.object.l2.HC.RData")
load("Cellchat.new.Seurat.object.l2.PC.RData")
load("Cellchat.new.Seurat.object.l2.WC.RData")
load("cellchat.new.l2.HC.DCE.RData")
load("cellchat.new.l2.HC.DRE.RData")
load("cellchat.new.l2.DCE.DRE.RData")
cellchatpathways <- unique(c(cellchat0@netP$pathways, cellchat1@netP$pathways,
                             cellchat2@netP$pathways, cellchatl2hc@netP$HC$pathways,
                             cellchatl2hc@netP$DCE$pathways, cellchatl2hr@netP$HC$pathways,
                             cellchatl2hr@netP$DRE$pathways, cellchatl2cr@netP$DCE$pathways,
                             cellchatl2cr@netP$DRE$pathways))
cellchatcompositelist <- c(cellchatl2hc, cellchatl2hr, cellchatl2cr)
names(cellchatcompositelist) <- c("HC vs DCE", "HC vs DRE", "DCE vs DRE")
cell_meta$TCR.cloneSize <- newcloneType

#Set factor levels of input objects and variables
l1_levels <- c("B", "Mono", "CD4 T", "CD8 T", "NK", "DC", "other T", "other")
l2_levels <- c("B naive", "B intermediate", "B memory", "Plasmablast", "CD14 Mono", "CD16 Mono",
               "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "Treg", "CD8 Naive",
               "CD8 TCM", "CD8 TEM", "gdT", "MAIT", "dnT", "CD4 Proliferating",
               "CD8 Proliferating", "NK", "NK_CD56bright", "NK Proliferating",
               "cDC2", "pDC", "Platelet", "HSPC", "ILC", "ASDC")
cell_meta$predicted.celltype.l1 <- factor(cell_meta$predicted.celltype.l1,
                                          levels = l1_levels)
cell_meta$predicted.celltype.l2 <- factor(cell_meta$predicted.celltype.l2,
                                          levels = l2_levels)

#Edit MC/WC to DCE and MR/PC to DRE in metadata
cell_meta$Condition <- str_replace(cell_meta$Condition, "MC", "DCE")
cell_meta$Condition <- str_replace(cell_meta$Condition, "MR", "DRE")
l1_l2$Condition <- str_replace(l1_l2$Condition, "MC", "DCE")
l1_l2$Condition <- str_replace(l1_l2$Condition, "MR", "DRE")
lis <- list(c("HC", "DCE"), c("HC", "DRE"), c("DCE", "DRE"))
MASCinput$Group <- str_replace(MASCinput$Group, "MC", "DCE")
MASCinput$Group <- str_replace(MASCinput$Group, "MR", "DRE")
pair_df$exclude <- str_replace(pair_df$exclude, "MC", "DCE")
pair_df$exclude <- str_replace(pair_df$exclude, "MR", "DRE")
pair_df$exclude <- factor(pair_df$exclude, levels = c("HC vs DCE", "HC vs DRE", "DRE vs DCE"))
PCWCL1$Group <- str_replace(PCWCL1$Group, "MC", "DCE")
PCWCL1$Group <- str_replace(PCWCL1$Group, "MR", "DRE")
l1_l2$Condition <- str_replace(l1_l2$Condition, "MC", "DCE")
l1_l2$Condition <- str_replace(l1_l2$Condition, "MR", "DRE")
summarized_data$exclude <- str_replace(summarized_data$exclude, "MC", "DCE")
summarized_data$exclude <- str_replace(summarized_data$exclude, "MR", "DRE")
summarized_data$exclude <- factor(summarized_data$exclude,
                                  levels = c("HC vs DCE", "HC vs DRE", "DRE vs DCE"))