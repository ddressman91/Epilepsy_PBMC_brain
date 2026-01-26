library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(colorspace)
library(readxl)
library(janitor)

#Load Seurat object
setwd("~")
load("Seurat.object.with.TCR.antigen.specificity.RData")

#UMAP projection of L2 cell types
seurat@meta.data$l2.celltypes.clean <- case_when(
  seurat@meta.data$predicted.celltype.l2 == "CD8Naive" ~ "CD8 naive",
  seurat@meta.data$predicted.celltype.l2 == "CD8TEM" ~ "CD8 TEM",
  seurat@meta.data$predicted.celltype.l2 == "Bnaive" ~ "B naive",
  seurat@meta.data$predicted.celltype.l2 == "CD14Mono" ~ "CD14 mono",
  seurat@meta.data$predicted.celltype.l2 == "cDC2" ~ "cDC2",
  seurat@meta.data$predicted.celltype.l2 == "CD16Mono" ~ "CD16 mono",
  seurat@meta.data$predicted.celltype.l2 == "Bmemory" ~ "B memory",
  seurat@meta.data$predicted.celltype.l2 == "CD4TEM" ~ "CD4 TEM",
  seurat@meta.data$predicted.celltype.l2 == "NK" ~ "NK CD56dim",
  seurat@meta.data$predicted.celltype.l2 == "CD4TCM" ~ "CD4 TCM",
  seurat@meta.data$predicted.celltype.l2 == "gdT" ~ "gdT",
  seurat@meta.data$predicted.celltype.l2 == "CD4Naive" ~ "CD4 naive",
  seurat@meta.data$predicted.celltype.l2 == "CD8TCM" ~ "CD8 TCM",
  seurat@meta.data$predicted.celltype.l2 == "Bintermediate" ~ "B intermediate",
  seurat@meta.data$predicted.celltype.l2 == "MAIT" ~ "MAIT",
  seurat@meta.data$predicted.celltype.l2 == "pDC" ~ "pDC",
  seurat@meta.data$predicted.celltype.l2 == "dnT" ~ "dnT",
  seurat@meta.data$predicted.celltype.l2 == "Treg" ~ "Treg",
  seurat@meta.data$predicted.celltype.l2 == "ILC" ~ "ILC",
  seurat@meta.data$predicted.celltype.l2 == "NKCD56bright" ~ "NK CD56bright",
  seurat@meta.data$predicted.celltype.l2 == "CD4CTL" ~ "CD4 CTL",
  seurat@meta.data$predicted.celltype.l2 == "NKProliferating" ~ "NK proliferating",
  seurat@meta.data$predicted.celltype.l2 == "HSPC" ~ "HSPC",
  seurat@meta.data$predicted.celltype.l2 == "ASDC" ~ "ASDC",
  seurat@meta.data$predicted.celltype.l2 == "CD4Proliferating" ~ "CD4 proliferating",
  seurat@meta.data$predicted.celltype.l2 == "CD8Proliferating" ~ "CD8 proliferating",
  seurat@meta.data$predicted.celltype.l2 == "Plasmablast" ~ "Plasmablast",
  seurat@meta.data$predicted.celltype.l2 == "Platelet" ~ "Platelet")
Idents(seurat) <- "l2.celltypes.clean"
DimPlot(seurat, label = T, repel = T) + theme(legend.position = "none") +
  ggtitle("UMAP projection of immune cell types")

#Dot plot of cell type lineage marker expression
seurat@meta.data$l2.celltypes.clean <- factor(seurat@meta.data$l2.celltypes.clean,
                                                 levels = rev(c("B naive", "B intermediate", "B memory",
                                                            "Plasmablast", "CD4 naive", "CD4 TCM",
                                                            "CD4 TEM", "CD4 CTL", "Treg", "CD8 naive", "CD8 TCM",
                                                            "CD8 TEM", "CD14 mono", "CD16 mono", "cDC2",
                                                            "pDC", "gdT", "MAIT", "NK CD56dim", "NK CD56bright",
                                                            "NK proliferating", "HSPC", "Platelet", "dnT",
                                                            "ILC", "CD4 proliferating", "CD8 proliferating", "ASDC")))
Idents(seurat) <- "l2.celltypes.clean"
DotPlot(seurat, features = c("CD3E", "CD4", "IL7R", "CCR7", "CD8A", "NOSIP", "GZMA", "NKG7",
                            "GNLY", "FOXP3", "TRAV1-2", "TRGC1", "TRGV9", "TRDV2", "CD79A",
                            "MS4A1", "POU2AF1", "HLA-DRB1", "CST3", "FCER1A", "LYZ", "S100A9",
                            "CD14", "FCGR3A", "FCER1G"),
        idents = c("B naive", "B intermediate", "B memory",
                   "Plasmablast", "CD4 naive", "CD4 TCM",
                   "CD4 TEM", "CD4 CTL", "Treg", "CD8 naive", "CD8 TCM",
                   "CD8 TEM", "CD14 mono", "CD16 mono", "cDC2",
                   "pDC", "gdT", "MAIT", "NK CD56dim", "NK CD56bright",
                   "NK proliferating", "HSPC", "Platelet"),
        dot.min = 0.1, col.min = 0, assay = "SCT") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()) +
  ggtitle("Dot plot of selected cell type markers")

#Donut plot of cell type proportions for L1 and L2 cell types, excluding
## cell types with <100 cells
celltypes <- read_xlsx("Supplementary table 1.xlsx", sheet = 2)
celltypes <- celltypes[nrow(celltypes),]
celltypes <- as.data.frame(t(celltypes))
celltypes <- row_to_names(celltypes, row_number = 1)
celltypes$Cell.Type <- rownames(celltypes)
celltypes <- celltypes[1:(nrow(celltypes)-1),]
colnames(celltypes)[1] <- "Count"
celltypes$Count <- as.numeric(celltypes$Count)
celltypes$Proportion <- celltypes$Count / sum(celltypes$Count)
celltypes$Cell.Type <- as.character(celltypes$Cell.Type)
celltypes$Percentage <- as.numeric(format(round(celltypes$Proportion * 100, 2), nsmall = 2))
excludecelltypes <- c("Plasmablast", "dnT", "HSPC", "Platelet",
                      "ASDC", "CD4 Proliferating", "CD8 Proliferating", "cDC1", "ILC",
                      "NK Proliferating", "NK_CD56bright")
celltypes$Cell.Type[celltypes$Cell.Type %in% excludecelltypes] <- "Other"
celltypes$Cell.Type <- factor(celltypes$Cell.Type,
                              levels = c("B naive", "B intermediate", "B memory",
                                         "CD4 Naive", "CD4 TCM", "CD4 TEM",
                                         "CD4 CTL", "Treg", "CD8 Naive", "CD8 TCM",
                                         "CD8 TEM", "cDC2", "pDC", "CD14 Mono", "CD16 Mono",
                                         "NK", "gdT", "MAIT", "Other"))

celltypes$Percentage[celltypes$Cell.Type == "Other"] <-
  sum((celltypes %>% filter(Cell.Type == "Other"))$Percentage)
celltypes$Proportion[celltypes$Cell.Type == "Other"] <-
  sum((celltypes %>% filter(Cell.Type == "Other"))$Proportion)
celltypes$Count[celltypes$Cell.Type == "Other"] <-
  sum((celltypes %>% filter(Cell.Type == "Other"))$Count)
celltypes <- unique(celltypes)
colors <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A",
  "#FF7F00", "gold1", "skyblue2", "#FB9A99", "palegreen2",
  "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1",
  "deeppink1", "blue1", "steelblue4", "darkturquoise")
celltypes$labels <- paste(celltypes$Cell.Type, " (", celltypes$Percentage, "%)", sep = "")
rownames(celltypes) <- c(1:nrow(celltypes))
labels <- rev(celltypes$labels)

ggplot(celltypes, aes(y = Proportion, fill = Cell.Type, x = 2)) +
  geom_col() + coord_polar("y") + theme_void() + xlim(c(0.5,2.5)) +
  annotate(geom = "text", x = 0.5, y = 0, label = "84K", size = 16, color = "black") +
  scale_fill_manual(values = colors, name = "Cell type", labels = labels) +
  guides(fill = guide_legend(ncol = 2))

#MASC plot of DCE and DRE vs HC
mascdce <- read_excel("Supplementary table 2.xlsx", sheet = 1)
mascdre <- read_excel("Supplementary table 2.xlsx", sheet = 2)
mascdce$Group <- "DCE"
mascdre$Group <- "DRE"
mascdce$cluster[mascdce$cluster == "NK"] <- "NKCD56dim"
mascdre$cluster[mascdre$cluster == "NK"] <- "NKCD56dim"
MASCinput <- bind_rows(mascdce, mascdre)
MASCinput$CI.upper <- as.numeric(MASCinput$CI.upper)
MASCinput <- MASCinput %>% filter(size > 100, CI.upper < 45, cluster != "HSPC")
MASCinput$Group <- factor(MASCinput$Group, levels = c("DCE", "DRE"))
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
           ymax = log2(CI.upper), shape = Group, col = Significance)) +
  geom_pointrange(position = position_dodge2(width = 0.7)) +
  xlab("Cell type") + ylab("log2(Odds Ratio)") + #ylim(c(-8,5.1)) +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("red", "black", "blue"))
