library(dplyr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scales)
library(readxl)

#Read in files with DGE output
setwd("~")
sheet_names <- excel_sheets("Supplementary table 3.xlsx")
DCEvsHC <- map(sheet_names,
               ~ read_excel("Supplementary table 3.xlsx", sheet = .x)) %>% 
  set_names(sheet_names)
sheet_names <- excel_sheets("Supplementary table 4.xlsx")
DREvsHC <- map(sheet_names,
               ~ read_excel("Supplementary table 4.xlsx", sheet = .x)) %>% 
  set_names(sheet_names)

#Run GSEA for DCE vs HC
DCEvsHCgseainput <- vector("list", 25)
names(DCEvsHCgseainput) <- names(DCEvsHC)
DCEvsHCgseaoutput <- vector("list", 25)
names(DCEvsHCgseaoutput) <- names(DCEvsHC)
for (i in 1:length(DCEvsHC)) {
  DCEvsHCgseainput[[i]] = DCEvsHC[[i]][,"logFC"]
  names(DCEvsHCgseainput[[i]]) = as.character(DCEvsHC[[i]][,"gene"])
  DCEvsHCgseainput[[i]] = sort(DCEvsHCgseainput[[i]], decreasing = TRUE)
  DCEvsHCgseaoutput[[i]] <- gseGO(geneList = DCEvsHCgseainput[[i]],
                             ont = "BP", OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",
                             exponent = 1, minGSSize = 15, maxGSSize = 200, eps = 1e-10,
                             pvalueCutoff = 0.05, pAdjustMethod = "fdr", by = "fgsea")
  DCEvsHCgseaoutput[[i]] <- DCEvsHCgseaoutput[[i]]@result
  if (nrow(DCEvsHCgseaoutput[[i]])>0) {
    DCEvsHCgseaoutput[[i]]$CellType <- names(DCEvsHCgseaoutput)[i]
      }
}

save(DCEvsHCgseaoutput, file = "GSEA.PBMCs.DCE.vs.HC.RData")

#Run GSEA for DRE vs HC
DREvsHCgseainput <- vector("list", 25)
names(DREvsHCgseainput) <- names(DREvsHC)
DREvsHCgseaoutput <- vector("list", 25)
names(DREvsHCgseaoutput) <- names(DREvsHC)
for (i in 1:length(DREvsHC)) {
  DREvsHCgseainput[[i]] = DREvsHC[[i]][,"logFC"]
  names(DREvsHCgseainput[[i]]) = as.character(DREvsHC[[i]][,"gene"])
  DREvsHCgseainput[[i]] = sort(DREvsHCgseainput[[i]], decreasing = TRUE)
  DREvsHCgseaoutput[[i]] <- gseGO(geneList = DREvsHCgseainput[[i]],
                                  ont = "BP", OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",
                                  exponent = 1, minGSSize = 15, maxGSSize = 200, eps = 1e-10,
                                  pvalueCutoff = 0.05, pAdjustMethod = "fdr", by = "fgsea")
  DREvsHCgseaoutput[[i]] <- DREvsHCgseaoutput[[i]]@result
  if (nrow(DREvsHCgseaoutput[[i]])>0) {
    DREvsHCgseaoutput[[i]]$CellType <- names(DREvsHCgseaoutput)[i]
  }
}

save(DREvsHCgseaoutput, file = "GSEA.PBMCs.DRE.vs.HC.RData")