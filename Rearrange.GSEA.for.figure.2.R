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
                                           "NK","other"))

combinedpathways <- unique(gseacombined$Description)
combinedpathways <- combinedpathways[str_detect(combinedpathways, "regulation") == F]
combinedpathways <- combinedpathways[c(1,3,69,75,77,86,99,104,115,122,126,131,137,141,144,
                                       145,146,150,156,162,176,181,183,189,191,226,230,231,
                                       244,300,323,326,356,357,360,376,383,389,390,404,436,
                                       477,499,501,502,526,529,531,536,606,643,655,691,695)]
combinedpathways <- rev(sort(combinedpathways))
combinedconcise <- combinedpathways[c(1,2,3,4,5,8,9,10,12,13,15,16,17,19,22,24,25,27,28,
                                      30,31,34,35,36,37,38,41,42,44,45,46,49,50,51,53,54)]

save(gseacombined, file = "GSEA.pathways.compiled.for.plotting.RData")
save(combinedpathways, file = "GSEA.immune.pathways.for.plotting.RData")
save(combinedconcise, file = "GSEA.immune.pathways.concise.for.plotting.RData")

#Brain GSEA plot
load("C:/Users/dalli/Documents/Epilepsy data/Brain snRNAseq/Braincells.GSEA.epilepsy.vs.control.04282025.RData")

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
       aes(x = CellType, y = Description, size = -log10(adj.P.Val), col = NES)) +
  geom_point() + theme_bw() + ggtitle("Selected GSEA pathways in TLE vs control") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(plot.title = element_text(hjust = 0.75)) +
  theme(axis.title.y = element_blank()) +
  guides(size = guide_legend(title = "-log10(p)")) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,1.1),"cm")) +
  scale_y_discrete(labels = wrap_format(70))