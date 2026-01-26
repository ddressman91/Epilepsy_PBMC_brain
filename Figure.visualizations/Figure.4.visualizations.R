library(dplyr)
library(stringr)
library(EnhancedVolcano)
library(scales)
library(readxl)

setwd("~")
dcedpe <- read_excel("Supplementary table 10.xlsx", sheet = 2)
dredpe <- read_excel("Supplementary table 10.xlsx", sheet = 3)
gseadcevshc <- read_excel("Supplementary table 10.xlsx", sheet = 5)
gseadrevshc <- read_excel("Supplementary table 10.xlsx", sheet = 6)

#GSEA plot comparing DRE vs HC to DCE vs HC
allgsea <- bind_rows(gseadcevshc %>% mutate(Comparison = "DCE vs HC"),
                     olinkgseadre %>% mutate(Comparison = "DRE vs HC"))
allgsea <- allgsea[c(1,3,8,13,16,21,23,29,33,36,44,46,51),]

ggplot(allgsea, aes(x = Comparison, y = Description, size = -log10(p.adjust), col = NES)) +
  geom_point() + theme_bw() + ylab(label = " ") +
  ggtitle("Selected pathways in DRE and DCE vs HC") +
  scale_y_discrete(labels = wrap_format(40))+
  theme(plot.title = element_text(hjust = 0.65)) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(axis.title.x = element_blank()) +
  scale_color_gradient(low = "#3F7CAC", high = "#56B1F7")

#Volcano plots of DCE vs HC and DRE vs HC results
dcelabel <- (dcedpe %>% filter(Adjusted_pval<0.05))$Assay
drelabel <- (dredpe %>% filter(Adjusted_pval<0.05))$Assay

EnhancedVolcano(dcedpe, lab = dcedpe$Assay,
                x = "estimate", y = "Adjusted_pval",
                xlab = "Effect size (DCE vs HC)", selectLab = dcelabel,
                title = "Proteomic changes in DCE patients", xlim = c(-3,3),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, dreutoff = 0.05,
                legendPosition = "None", caption = NULL, ylim = c(0,2),
                drawConnectors = T, arrowheads = F, max.overlaps = 75) +
  theme(plot.title = element_text(face = "plain"))

EnhancedVolcano(dredpe, lab = dredpe$Assay,
                x = "estimate", y = "Adjusted_pval", xlim = c(-3.5,2.5),
                xlab = "Effect size (DRE vs HC)", selectLab = drelabel,
                title = "Proteomic changes in DRE patients", ylim = c(0,1.5),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, dreutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 75) +
  theme(plot.title = element_text(face = "plain"))
