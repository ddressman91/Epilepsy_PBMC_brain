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
library(ggplot2)
library(stringr)
library(CellChat)
library(cowplot)
library(scRepertoire)
library(scales)
library(ggiraph)

load("Alternate.DimPlot.function.RData")
load("Alternate.FeaturePlot.function.RData")
load("Alternate.VlnPlot.function.RData")
load("TLEbrain.Seurat.metadata.downsampled.RData")
load("TLEbrain.Seurat.RNA.data.matrix.downsampled.max50plus.RData")
load("TLEbrain.Seurat.UMAP.coordinates.downsampled.RData")
mgclust <- read.csv("Microglia.cluster.markers.epilepsy.brain.all.samples.csv")
colnames(mgclust)[1] <- "Gene"
mgclust$Gene <- gsub("\\.*", "", mgclust$Gene)
load("Braincells.DGE.epilepsy.vs.control.04282025.RData")
for (i in 1:length(braincellsdds)) {
  braincellsdds[[i]]$gene <- rownames(braincellsdds[[i]])
}
load("Brain.CellChat.epilepsy.downsampled.control.merged.RData")
load("Brain.CellChat.epilepsy.downsampled.01082026.RData")
load("Brain.CellChat.ctrl.04232025.RData")
brainpathways <- unique(c(cellchatctrl@netP$pathways, cellchatepilepsy@netP$pathways,
                          cellchatepctrl@netP$Control$pathways,
                          cellchatepctrl@netP$Epilepsy$pathways))
cell_meta$CellType <- str_replace(cell_meta$CellType, "CD8\\+ T", "CD8T")
cell_meta$CellType <- str_replace(cell_meta$CellType, "CD4\\+ T", "CD4T")

braincelltypes <- unique(cell_meta$CellType)
braincelltypes <- factor(braincelltypes, levels = c("MG1", "MG2", "MG3", "MG4",
                                                    "MG5", "MG6", "MG7", "MG8",
                                                    "MG9", "Monocytes", "CD4T",
                                                    "CD8T", "gdT", "NK", "B cells", "Other"))

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
