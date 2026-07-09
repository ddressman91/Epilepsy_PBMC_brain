library(dplyr)
library(stringr)
library(ggplot2)
library(EnhancedVolcano)

#Load in GSEA and DGE files
setwd("~")
sheet_names <- excel_sheets("Supplementary table 5.xlsx")
dcegsealist <- map(sheet_names,
                ~ read_excel("Supplementary table 5.xlsx", sheet = .x)) %>% 
  set_names(sheet_names)

sheet_names <- excel_sheets("Supplementary table 6.xlsx")
dregsealist <- map(sheet_names,
                   ~ read_excel("Supplementary table 6.xlsx", sheet = .x)) %>% 
  set_names(sheet_names)

sheet_names <- excel_sheets("Supplementary table 3.xlsx")
dcedgelist <- map(sheet_names,
                   ~ read_excel("Supplementary table 3.xlsx", sheet = .x)) %>% 
  set_names(sheet_names)

sheet_names <- excel_sheets("Supplementary table 4.xlsx")
dredgelist <- map(sheet_names,
                   ~ read_excel("Supplementary table 4.xlsx", sheet = .x)) %>% 
  set_names(sheet_names)

l1celltypes <- names(dcegsealist)[c(1:7)]
l2celltypes <- names(dcegsealist)[c(8:25)]

#Examine convergent genes and pathways in DCE and DRE comparisons to HC, or
#genes and pathways that are significant with the same direction of effect
#in the same cell type in both DCE vs HC and DRE vs HC
dgeconvup <- vector("list", 25)
names(dgeconvup) <- names(dcegsealist)
dgeconvdown <- vector("list", 25)
names(dgeconvdown) <- names(dcegsealist)
gseaconvup <- vector("list", 25)
names(gseaconvup) <- names(dcegsealist)
gseaconvdown <- vector("list", 25)
names(gseaconvdown) <- names(dcegsealist)
for (i in 1:25) {
  dgeconvup[[i]] <- intersect((dcedgelist[[i]]%>%filter(adj.P.Val<0.05,logFC>0))$gene,
                                (dredgelist[[i]]%>%filter(adj.P.Val<0.05,logFC>0))$gene)
  dgeconvdown[[i]] <- intersect((dcedgelist[[i]]%>%filter(adj.P.Val<0.05,logFC<0))$gene,
                                (dredgelist[[i]]%>%filter(adj.P.Val<0.05,logFC<0))$gene)
  gseaconvup[[i]] <- intersect((dcegsealist[[i]] %>% filter(NES>0))$Description,
                               (dregsealist[[i]] %>% filter(NES>0))$Description)
  gseaconvdown[[i]] <- intersect((dcegsealist[[i]] %>% filter(NES<0))$Description,
                               (dregsealist[[i]] %>% filter(NES<0))$Description)
}

#Examine specific genes and pathways found only in DCE vs HC or DRE vs HC
drespecific <- vector("list", 25)
dcespecific <- vector("list", 25)
names(drespecific) <- names(dcegsealist)
names(dcespecific) <- names(dcegsealist)
for (i in 1:25) {
  dregenes <- setdiff((dredgelist[[i]]%>%filter(adj.P.Val<0.05))$gene,(dcedgelist[[i]]%>%filter(adj.P.Val<0.05))$gene)
  dcegenes <- setdiff((dcedgelist[[i]]%>%filter(adj.P.Val<0.05))$gene,(dredgelist[[i]]%>%filter(adj.P.Val<0.05))$gene)
  drespecific[[i]] <- dredgelist[[i]] %>% filter(gene %in% dregenes)
  dcespecific[[i]] <- dcedgelist[[i]] %>% filter(gene %in% dcegenes)
}

dregseaspecific <- vector("list", 25)
dcegseaspecific <- vector("list", 25)
names(dregseaspecific) <- names(dcegsealist)
names(dcegseaspecific) <- names(dcegsealist)
for (i in 1:25) {
  drepathways <- setdiff((dregsealist[[i]])$Description,(dcegsealist[[i]])$Description)
  dcepathways <- setdiff((dcegsealist[[i]])$Description,(dregsealist[[i]])$Description)
  dregseaspecific[[i]] <- dregsealist[[i]] %>% filter(Description %in% drepathways)
  dcegseaspecific[[i]] <- dcegsealist[[i]] %>% filter(Description %in% dcepathways)
}

alldcespecific <- c()
for (i in 1:length(dcegseaspecific)) {
  alldcespecific <- rbind(alldcespecific, dcegseaspecific[[i]])
}

alldrespecific <- c()
for (i in 1:length(dregseaspecific)) {
  alldrespecific <- rbind(alldrespecific, dregseaspecific[[i]])
}

#Examine genes and pathways significant in both DCE vs HC and DRE vs HC in
#the same cell type, regardless of direction of effect
alldcegseaconvup <- c()
for (i in 1:25) {
  if (nrow(dcegsealist[[i]])>0) {
    alldcegseaconvup <- bind_rows(alldcegseaconvup, dcegsealist[[i]] %>%
                                    filter(Description %in% gseaconvup[[i]]))
  }
}
alldcegseaconvdown <- c()
for (i in 1:25) {
  if (nrow(dcegsealist[[i]])>0) {
    alldcegseaconvdown <- bind_rows(alldcegseaconvdown, dcegsealist[[i]] %>%
                                    filter(Description %in% gseaconvdown[[i]]))
  }
}
alldregseaconvup <- c()
for (i in 1:25) {
  if (nrow(dregsealist[[i]])>0) {
    alldregseaconvup <- bind_rows(alldregseaconvup, dregsealist[[i]] %>%
                                    filter(Description %in% gseaconvup[[i]]))
  }
}
alldregseaconvdown <- c()
for (i in 1:25) {
  if (nrow(dregsealist[[i]])>0) {
    alldregseaconvdown <- bind_rows(alldregseaconvdown, dregsealist[[i]] %>%
                                    filter(Description %in% gseaconvdown[[i]]))
  }
}

gseaconvplot <- bind_rows(alldcegseaconvup %>% mutate(Group = "DCE", Direction = "Up in epilepsy"),
                          alldcegseaconvdown %>% mutate(Group = "DCE", Direction = "Down in epilepsy"),
                          alldregseaconvup %>% mutate(Group = "DRE", Direction = "Up in epilepsy"),
                          alldregseaconvdown %>% mutate(Group = "DRE", Direction = "Down in epilepsy"))

#Make a composite GSEA table
gseaconvplot$Group <- "All epilepsy"
alldcespecific$Group <- "DCE only"
alldrespecific$Group <- "DRE only"
alldcespecific$Direction <- case_when(alldcespecific$NES > 0 ~ "Up in epilepsy",
                                     alldcespecific$NES < 0 ~ "Down in epilepsy")
alldrespecific$Direction <- case_when(alldrespecific$NES > 0 ~ "Up in epilepsy",
                                     alldrespecific$NES < 0 ~ "Down in epilepsy")

gseacombined <- bind_rows(gseaconvplot, alldcespecific, alldrespecific)
gseacombined$Direction <- factor(gseacombined$Direction, levels = c("Up in epilepsy",
                                                                    "Down in epilepsy"))
gseacombined$Group <- factor(gseacombined$Group, levels = c("DCE only", "DRE only", "All epilepsy"))
gseacombined$CellType <- factor(gseacombined$CellType,
                                levels = c("B","Bnaive","Bintermediate","Bmemory",
                                           "CD4T","CD4Naive","CD4TCM","CD4TEM","CD4CTL","Treg",
                                           "CD8T","CD8Naive","CD8TCM","CD8TEM","otherT","gdT","MAIT",
                                           "Mono","CD14Mono","CD16Mono","DC","cDC2","pDC",
                                           "NK","NKCD56dim", "NKCD56bright", "other"))

combinedpathways <- unique(gseacombined$Description)
combinedpathways <- combinedpathways[str_detect(combinedpathways, "regulation") == F]
combinedconcise <- combinedpathways[c(3,4,12,15,18,19,21,23,25,33:35,37,41,42,45,62,
                                      79,80,89,91,101,133,135,229:231,247,259,260,
                                      267,276,285)]

save(gseacombined, file = "GSEA.pathways.compiled.for.plotting.RData")
save(combinedconcise, file = "GSEA.immune.pathways.concise.for.plotting.RData")
