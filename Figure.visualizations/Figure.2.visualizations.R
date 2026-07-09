library(dplyr)
library(stringr)
library(ggplot2)
library(EnhancedVolcano)
library(scales)
library(readxl)
library(rlang)
library(purrr)

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
ggplot(gseacombined %>% filter(Description %in% combinedconcise),
       aes(x = CellType, y = Description, shape = Direction, col = Group, size = abs(NES))) +
  geom_point() + theme_bw() + ggtitle("Shared and unique pathways by treatment response") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(plot.title = element_text(hjust = 0.75)) +
  theme(axis.title.y = element_blank()) +
  guides(size = guide_legend(title = "abs(NES)")) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,1.1),"cm")) +
  scale_y_discrete(labels = wrap_format(70)) +
  scale_shape_manual(values = c("\u25B2","\u25BC")) +
  scale_color_manual(values = c("#FFA500","#7EC0EE","black"))

write.csv(gseacombined %>% filter(Description %in% combinedconcise),
          file = "GSEA.table.for.figure.2.plot.csv", row.names = F)

#Make volcano plots of selected cell types
blabel <- (dcedgelist$B %>% filter(adj.P.Val < 0.05))$gene
blabel <- blabel[str_detect(blabel, "AC[0-9]|AL[0-9]|HIST|RPS|RPL|MT-") == F]
EnhancedVolcano(dcedgelist$B, lab = dcedgelist$B$gene,
                selectLab = blabel,
                x = "logFC", y = "adj.P.Val", xlim = c(-4,3),
                xlab = "Log2FC (DCE vs HC)", ylab = "-Log10(p.adjusted)",
                title = "B cell DEGs in DCE patients", ylim = c(0,8),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 30) +
  theme(plot.title = element_text(face = "plain"))

cd8tlabel <- (dredgelist$CD8T %>% filter(adj.P.Val < 0.05))$gene
cd8tlabel <- cd8tlabel[str_detect(cd8tlabel,
                                          "AC[0-9]|AL[0-9]|HIST|RPS|RPL|MT-") == F]
EnhancedVolcano(dredgelist$CD8T, lab = dredgelist$CD8T$gene,
                selectLab = cd8tlabel,
                x = "logFC", y = "adj.P.Val", xlim = c(-3.5,3.5),
                xlab = "Log2FC (DRE vs HC)", ylab = "-Log10(p.adjusted)",
                title = "CD8+ T cell DEGs in DRE patients", ylim = c(0,5),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 30) +
  theme(plot.title = element_text(face = "plain"))

cd14monolabel <- (dredgelist$CD14Mono %>% filter(adj.P.Val < 0.05))$gene
cd14monolabel <- cd14monolabel[str_detect(cd14monolabel,
                                  "AC[0-9]|AL[0-9]|HIST|RPS|RPL|MT-") == F]
EnhancedVolcano(dredgelist$CD14Mono, lab = dredgelist$CD14Mono$gene,
                selectLab = cd14monolabel,
                x = "logFC", y = "adj.P.Val", xlim = c(-6,2.5),
                xlab = "Log2FC (DRE vs HC)", ylab = "-Log10(p.adjusted)",
                title = "CD14+ monocyte DEGs in DRE patients", ylim = c(0,6.5),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 30) +
  theme(plot.title = element_text(face = "plain"))

cd14monolabel <- (dcedgelist$CD14Mono %>% filter(adj.P.Val < 0.05))$gene
cd14monolabel <- cd14monolabel[str_detect(cd14monolabel,
                                          "AC[0-9]|AL[0-9]|HIST|RPS|RPL|MT-") == F]
EnhancedVolcano(dcedgelist$CD14Mono, lab = dcedgelist$CD14Mono$gene,
                selectLab = cd14monolabel,
                x = "logFC", y = "adj.P.Val", xlim = c(-6,2.5),
                xlab = "Log2FC (DCE vs HC)", ylab = "-Log10(p.adjusted)",
                title = "CD14+ monocyte DEGs in DCE patients", ylim = c(0,6),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 30) +
  theme(plot.title = element_text(face = "plain"))

#Volcano plots of DCE vs HC and DRE vs HC results from proteomics
dcedpe <- read_excel("Supplementary table 13.xlsx", sheet = 2)
dcedpeposthoc <- read_excel("Supplementary table 13.xlsx", sheet = 3)
dredpe <- read_excel("Supplementary table 13.xlsx", sheet = 4)
dredpeposthoc <- read_excel("Supplementary table 13.xlsx", sheet = 5)

dcelabel <- (dcedpe %>% filter(Adjusted_pval<0.05))$Assay
drelabel <- (dredpe %>% filter(Adjusted_pval<0.05))$Assay

dcedpe <- left_join(dcedpe, dcedpeposthoc %>% dplyr::select(Assay,estimate),
                    by = "Assay")
dredpe <- left_join(dredpe, dredpeposthoc %>% dplyr::select(Assay,estimate),
                    by = "Assay")

EnhancedVolcano(dcedpe, lab = dcedpe$Assay,
                x = "estimate", y = "Adjusted_pval",
                xlab = "Effect size (DCE vs HC)", selectLab = dcelabel,
                title = "Proteomic changes in DCE patients", xlim = c(-5.5,8),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL, ylim = c(0,8),
                drawConnectors = T, arrowheads = F, max.overlaps = 75) +
  theme(plot.title = element_text(face = "plain"))

EnhancedVolcano(dredpe, lab = dredpe$Assay,
                x = "estimate", y = "Adjusted_pval", xlim = c(-5,6),
                xlab = "Effect size (DRE vs HC)", selectLab = drelabel,
                title = "Proteomic changes in DRE patients", ylim = c(0,2.5),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 75) +
  theme(plot.title = element_text(face = "plain"))

#GSEA plot comparing DRE, DCE, and HC in proteomics
gseadcevshc <- read_excel("Supplementary table 13.xlsx", sheet = 8)
gseadrevshc <- read_excel("Supplementary table 13.xlsx", sheet = 9)
gseadrevsdce <- read_excel("Supplementary table 13.xlsx", sheet = 10)

allgsea <- bind_rows(gseadcevshc %>% mutate(Comparison = "DCE vs HC"),
                     gseadrevshc %>% mutate(Comparison = "DRE vs HC"),
                     gseadrevsdce %>% mutate(Comparison = "DRE vs DCE"))
pathways <- allgsea$Description[c(1:3,13,25,27:29,32,34,44,54,57,65,72:74,87)]
allgsea$Comparison <- factor(allgsea$Comparison,
                             levels = c("DCE vs HC", "DRE vs HC", "DRE vs DCE"))

ggplot(allgsea %>% filter(Description %in% pathways),
       aes(x = Comparison, y = Description, size = -log10(p.adjust), col = NES)) +
  geom_point() + theme_bw() + ylab(label = " ") +
  scale_y_discrete(labels = wrap_format(50)) +
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  theme(axis.title.x = element_blank())
