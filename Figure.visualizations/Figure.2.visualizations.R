library(dplyr)
library(stringr)
library(ggplot2)
library(EnhancedVolcano)
library(scales)

#Load in GSEA and DGE files
setwd("~")
sheet_names <- excel_sheets("Supplementary table 3.xlsx")
dcedgelist <- map(sheet_names,
                  ~ read_excel("Supplementary table 3.xlsx", sheet = .x)) %>% 
  set_names(sheet_names)

sheet_names <- excel_sheets("Supplementary table 4.xlsx")
dredgelist <- map(sheet_names,
                  ~ read_excel("Supplementary table 4.xlsx", sheet = .x)) %>% 
  set_names(sheet_names)

load("GSEA.pathways.compiled.for.plotting.RData") #See Rearrange.GSEA.for.figure.2.R
load("GSEA.immune.pathways.concise.for.plotting.RData") #See Rearrange.GSEA.for.figure.2.R

#Reproduce GSEA figure 2 plot
ggplot(gseacombined %>% filter(Description %in% combinedconcise,
                               !CellType %in% c("B", "CD8T", "otherT", "Mono", "DC", "NK")),
       aes(x = CellType, y = Description, shape = Direction, col = Group)) +
  geom_point(size = 3) + theme_bw() + ggtitle("Shared and unique pathways by treatment response") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(plot.title = element_text(hjust = 0.75)) +
  theme(axis.title.y = element_blank()) +
  guides(size = guide_legend(title = "-log10(p)")) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,1.1),"cm")) +
  scale_y_discrete(labels = wrap_format(70)) +
  scale_shape_manual(values = c("\u25B2","\u25BC")) +
  scale_color_manual(values = c("#FFA500","#7EC0EE","black"))

write.csv(gseacombined %>% filter(Description %in% combinedconcise),
          file = "GSEA.table.for.figure.2.plot.csv", row.names = F)

#Make volcano plots of selected cell types
cd14monolabel <- (dredgelist$CD14Mono %>% filter(adj.P.Val < 0.05))$gene
cd14monolabel <- cd14monolabel[str_detect(cd14monolabel,
                                  "AC[0-9]|AL[0-9]|HIST|RPS|RPL|MT-") == F]
EnhancedVolcano(dredgelist$CD14Mono, lab = dredgelist$CD14Mono$gene,
                selectLab = cd14monolabel,
                x = "logFC", y = "adj.P.Val", xlim = c(-5.2,2),
                xlab = "Log2FC (DRE vs HC)",
                title = "CD14+ mono DEGs in DRE patients", ylim = c(0,16),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 30) +
  theme(plot.title = element_text(face = "plain"))

cdc2label <- (dredgelist$cDC2 %>% filter(adj.P.Val < 0.05))$gene
cdc2label <- cdc2label[str_detect(cdc2label,
                                  "AC[0-9]|AL[0-9]|HIST|RPS|RPL|MT-") == F]
EnhancedVolcano(dredgelist$cDC2, lab = dredgelist$cDC2$gene,
                selectLab = cdc2label,
                x = "logFC", y = "adj.P.Val", xlim = c(-2.8,1.5),
                xlab = "Log2FC (DRE vs HC)",
                title = "cDC2 DEGs in DRE patients", ylim = c(0,10.1),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 50) +
  theme(plot.title = element_text(face = "plain"))

CD4Naivelist <- gseacombined %>% filter(Description %in% combinedconcise,
                                      CellType == "CD4Naive", Group == "DCE only")
CD4Naivecoregenes <- unique(unlist(str_split(CD4Naivelist$core_enrichment, pattern = "/")))
CD4Naivelabel <- (dcedgelist$CD4Naive %>% filter(adj.P.Val < 0.05))$gene
CD4Naivelabel <- CD4Naivelabel[str_detect(CD4Naivelabel,
                                      "AC[0-9]|AL[0-9]|HIST|RPS|RPL|MT-") == F]
CD4Naivelabel <- intersect(CD4Naivelabel, CD4Naivecoregenes)
EnhancedVolcano(dcedgelist$CD4Naive, lab = dcedgelist$CD4Naive$gene,
                selectLab = CD4Naivelabel,
                x = "logFC", y = "adj.P.Val", xlim = c(-2,2),
                xlab = "Log2FC (DCE vs HC)",
                title = "CD4+ Naive DEGs in DCE patients", ylim = c(0,8),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 75) +
  theme(plot.title = element_text(face = "plain"))

CD8Naivelist <- gseacombined %>% filter(Description %in% combinedconcise,
                                        CellType == "CD8Naive", Group == "DCE only")
CD8Naivecoregenes <- unique(unlist(str_split(CD8Naivelist$core_enrichment, pattern = "/")))
CD8Naivelabel <- (dcedgelist$CD8Naive %>% filter(adj.P.Val < 0.05))$gene
CD8Naivelabel <- CD8Naivelabel[str_detect(CD8Naivelabel,
                                          "AC[0-9]|AL[0-9]|HIST|RPS|RPL|MT-") == F]
CD8Naivelabel <- intersect(CD8Naivelabel, CD8Naivecoregenes)
EnhancedVolcano(dcedgelist$CD8Naive, lab = dcedgelist$CD8Naive$gene,
                selectLab = CD8Naivelabel,
                x = "logFC", y = "adj.P.Val", xlim = c(-2,2),
                xlab = "Log2FC (DCE vs HC)",
                title = "CD8+ Naive DEGs in DCE patients", ylim = c(0,6),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 75) +
  theme(plot.title = element_text(face = "plain"))

CD4TCMlabel <- (dcedgelist$CD4TCM %>% filter(adj.P.Val < 0.05))$gene
CD4TCMlabel <- CD4TCMlabel[str_detect(CD4TCMlabel,
                                      "AC[0-9]|AL[0-9]|HIST|RPS|RPL|MT-") == F]
EnhancedVolcano(dcedgelist$CD4TCM, lab = dcedgelist$CD4TCM$gene,
                selectLab = CD4TCMlabel,
                x = "logFC", y = "adj.P.Val", xlim = c(-2.5,2.2),
                xlab = "Log2FC (DCE vs HC)",
                title = "CD4+ TCM DEGs in DCE patients", ylim = c(0,6),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 50) +
  theme(plot.title = element_text(face = "plain"))

CD8TCMlabel <- (dcedgelist$CD8TCM %>% filter(adj.P.Val < 0.05))$gene
CD8TCMlabel <- CD8TCMlabel[str_detect(CD8TCMlabel,
                                "AC[0-9]|AL[0-9]|HIST|RPS|RPL|MT-") == F]
EnhancedVolcano(dcedgelist$CD8TCM, lab = dcedgelist$CD8TCM$gene,
                selectLab = CD8TCMlabel,
                x = "logFC", y = "adj.P.Val", xlim = c(-2,2),
                xlab = "Log2FC (DCE vs HC)",
                title = "CD8+ TCM DEGs in DCE patients", ylim = c(0,5.2),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 50) +
  theme(plot.title = element_text(face = "plain"))

gdTlabel <- (dcedgelist$gdT %>% filter(adj.P.Val < 0.05))$gene
gdTlabel <- gdTlabel[str_detect(gdTlabel,
                                  "AC[0-9]|AL[0-9]|HIST|RPS|RPL|MT-") == F]
EnhancedVolcano(dcedgelist$gdT, lab = dcedgelist$gdT$gene,
                selectLab = gdTlabel,
                x = "logFC", y = "adj.P.Val", xlim = c(-2.2,2),
                xlab = "Log2FC (DCE vs HC)", caption = NULL,
                title = "gdT DEGs in DCE patients", ylim = c(0,3.2), legendPosition = "None",
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendLabels = c("NS", "log2FC", "pval","pval + log2FC"),
                drawConnectors = T, arrowheads = F, max.overlaps = 75) +
  theme(plot.title = element_text(face = "plain"))

cdc2label <- (dcedgelist$cDC2 %>% filter(adj.P.Val < 0.05))$gene
cdc2label <- cdc2label[str_detect(cdc2label,
                                  "AC[0-9]|AL[0-9]|HIST|RPS|RPL|MT-") == F]
EnhancedVolcano(dcedgelist$cDC2, lab = dcedgelist$cDC2$gene,
                selectLab = cdc2label,
                x = "logFC", y = "adj.P.Val", xlim = c(-2,1.5),
                xlab = "Log2FC (DCE vs HC)",
                title = "cDC2 DEGs in DCE patients", ylim = c(0,5),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 50) +
  theme(plot.title = element_text(face = "plain"))
