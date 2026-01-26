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

#MASC plot
MASCoutput <- read_excel("Supplementary table 2.xlsx", sheet = 4)
MASCplot <- ggplot(MASCoutput, aes(x = CellType, y = OR, ymin = CI.lower,
                                   ymax = CI.upper, col = Status)) +
  geom_pointrange(position = position_dodge2(width = 0.7)) +
  xlab("Cell type") + ylab("Odds Ratio") + ylim(c(0,11)) +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  ggtitle("Immune cell enrichment in epilepsy brain") +
  scale_color_manual(values = c("red", "black", "blue"))

#CellChat
load("Brain.CellChat.epilepsy.control.merged.RData")
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

#Brain GSEA plot
sheet_names <- excel_sheets("Supplementary table 14.xlsx")
braingseaoutput <- map(sheet_names,
                ~ read_excel("Supplementary table 14.xlsx", sheet = .x)) %>% 
  set_names(sheet_names)
braingsea <- c()
for (i in c(1,3,4,6,9,11)) {
  braingsea <- bind_rows(braingsea,
                         braingseaoutput[[i]] %>% mutate(CellType = names(braingseaoutput)[i]))
}

brainpathways <- unique(braingsea$Description)
brainpathways <- brainpathways[str_detect(brainpathways,"regulation")==F]
brainpathways <- brainpathways[c(29,44,77,131,138,140,155,169,171,200,210,226,248,
                                 279,293,310,326,333,334,350,358,374,407,448,455,464,486,494,498,
                                 520,528,532,534,565,582,627,660,675,678,692,696,711,720,725)]
brainpathways <- brainpathways[c(1,2,4:20,22:30,32:36,38:44)] #Remove T cell intrinsic pathways

ggplot(braingsea %>% filter(Description %in% brainpathways),
       aes(x = CellType, y = Description, size = -log10(p.adjust), col = NES)) +
  geom_point() + theme_bw() + ggtitle("Selected GSEA pathways in TLE vs control") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(plot.title = element_text(hjust = 0.75)) +
  theme(axis.title.y = element_blank()) +
  guides(size = guide_legend(title = "-log10(p)")) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,1.1),"cm")) +
  scale_y_discrete(labels = wrap_format(70))

#Violin plot of GZMB in selected PBMC cell types
load("Seurat.object.with.TCR.antigen.specificity.RData")
Idents(seurat) <- "predicted.celltype.l2"
seurat@meta.data$condition <- case_when(seurat@meta.data$condition == 0 ~ "HC",
                                        seurat@meta.data$condition == 1 ~ "DCE",
                                        seurat@meta.data$condition == 2 ~ "DRE")
VlnPlot(seurat, features = "GZMB", idents = c("CD4CTL", "CD8TEM", "gdT", "MAIT", "NK", "Treg"),
        group.by = "predicted.celltype.l2", split.by = "condition")