library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(Seurat)
library(scRepertoire)
library(readxl)

#Load in data, UMAP plot of cell types
setwd("~")
load("Epilepsy.brain.celltypes.TCRs.04232025.RData")
Idents(tlebrain) <- "CellType"
DimPlot(tlebrain, label = T, repel = T)

#Dot plot of microglial clusters
microglia <- subset(tlebrain, subset = CellType %in% 
                      c("MG1", "MG2", "MG3", "MG4", "MG5", "MG6", "MG7", "MG8", "MG9"))
microglia <- SCTransform(microglia)
Idents(microglia) <- "CellType"
levels(Idents(microglia)) <- rev(c("MG1", "MG2", "MG3", "MG4", "MG5", "MG6", "MG7", "MG8", "MG9"))
figure5markers <- unique(c("IFIT1", "IFIT2", "IFIT3", "STAT1",
                           "DENND3", "PLCG2", "ATM",
                           "CD14", "CD83", "APOC1", "TSPO",
                           "GIMAP4", "UQCRB", "DAB2",
                           "CCL2", "CCL3", "CCL4", "IL1B", "TNF",
                           "CX3CR1", "P2RY12", "SELPLG", "SYK",
                           "CSF1R", "CSF2RA",
                           "FOS", "JUN", "IER2"))

DotPlot(microglia, features = figure5markers, assay = "SCT",
        idents = c("MG1", "MG2", "MG3", "MG4", "MG5", "MG6", "MG7", "MG8")) +
  theme(axis.text.x = element_text(angle=45,hjust=1))

#Brain GSEA plot
sheet_names <- excel_sheets("Supplementary table 22.xlsx")
braingseaoutput <- map(sheet_names,
                       ~ read_excel("Supplementary table 22.xlsx", sheet = .x)) %>% 
  set_names(sheet_names)
braingsea <- c()
for (i in 1:10) {
  if (nrow(braingseaoutput[[i]]) > 0) {
    braingsea <- bind_rows(braingsea,
                           braingseaoutput[[i]] %>% mutate(CellType = names(braingseaoutput)[i]))
  }
}

brainpathways <- unique(braingsea$Description)
brainpathways <- brainpathways[str_detect(brainpathways,"regulation")==F]
brainpathways <- brainpathways[c(33,39,108,113,125,187,202,242,247,251,270,312,328,
                                 336,337,345:347,349,351,354,361,366)]

ggplot(braingsea %>% filter(Description %in% brainpathways),
       aes(x = CellType, y = Description, size = -log10(p.adjust), col = NES)) +
  geom_point() + theme_bw() + ggtitle("Selected GSEA pathways in TLE vs control") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(plot.title = element_text(hjust = 0.75)) +
  theme(axis.title.y = element_blank()) +
  guides(size = guide_legend(title = "-log10(p)")) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,1.1),"cm")) +
  scale_y_discrete(labels = wrap_format(70))

#MASC plot
MASCoutput <- read_excel("Supplementary table 2.xlsx", sheet = 4)
MASCplot <- ggplot(MASCoutput %>% filter(Size > 600, CellType != "MG9"),
                   aes(x = CellType, y = log10(OR), ymin = log10(CI.lower),
                                   ymax = log10(CI.upper), col = Status)) +
  geom_pointrange(position = position_dodge2(width = 0.7)) +
  xlab("Cell type") + ylab("log10(Odds Ratio)") +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggtitle("Immune cell enrichment in epilepsy brain") +
  scale_color_manual(values = c("red", "black", "blue"))

#CellChat
load("Brain.CellChat.epilepsy.control.merged.RData")
load("Brain.CellChat.epilepsy.04232025.RData")
load("Brain.CellChat.ctrl.04232025.RData")
weight <- rankNet(cellchatepctrl, mode = "comparison", measure = "weight", sources.use = NULL,
                  targets.use = NULL, stacked = T, do.stat = TRUE,
                  title = "Cell-cell signaling pathways",
                  signaling = c("Prostaglandin", "SPP1", "ICAM", "GALECTIN", "CLEC",
                                "CXCL", "TGFb", "LAIR1", "CCL", "CSF", "IFN-II", "CX3C"))
pathways <- c("IFN-I", "NCAM", "CD86", "CD80", "TGFb", "IL2",
              "TNF", "IL4", "CCL", "TRAIL", "IL1", "SELL",
              "PECAM1", "SELPLG")
weight$theme$axis.text.y$colour <- c(rep("#00BFC4", 4), rep("#F8766D", 8))
weight

netVisual_aggregate(cellchatctrl, signaling = "CCL", layout = "circle",
                    remove.isolate = T)
netVisual_aggregate(cellchatepilepsy, signaling = "CCL", layout = "circle",
                    remove.isolate = T)

#Violin plot of GZMB in selected PBMC cell types
load("Seurat.object.with.TCR.antigen.specificity.RData")
Idents(seurat) <- "predicted.celltype.l2"
seurat@meta.data$condition <- case_when(seurat@meta.data$condition == 0 ~ "HC",
                                        seurat@meta.data$condition == 1 ~ "DCE",
                                        seurat@meta.data$condition == 2 ~ "DRE")
VlnPlot(seurat, features = "GZMB", idents = c("CD4CTL", "CD8TEM", "gdT", "MAIT", "NK", "Treg"),
        group.by = "predicted.celltype.l2", split.by = "condition",
        cols = c("#8B0000","#FFA500", "#7EC0EE"))
