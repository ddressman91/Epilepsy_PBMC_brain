library(Seurat)
library(CellChat)
library(patchwork)
library(future)
library(ComplexHeatmap)

#Load in Seurat object of PBMC data, subset by diagnosis
setwd("~")
load("Seurat.object.with.TCR.antigen.specificity.RData")
control <- subset(seurat, subset = condition == 0)
dce <- subset(seurat, subset = condition == 1)
dre <- subset(seurat, subset = condition == 2 & !PatientID %in% c(25,61))

#Create CellChat objects and add the database
cellchat0 <- createCellChat(object = control, group.by = "predicted.celltype.l2")
cellchat1 <- createCellChat(object = dce, group.by = "predicted.celltype.l2")
cellchat2 <- createCellChat(object = dre, group.by = "predicted.celltype.l2")
CellChatDB <- CellChatDB.human
cellchat0@DB <- CellChatDB
cellchat1@DB <- CellChatDB
cellchat2@DB <- CellChatDB

#Preprocessing of CellChat objects
cellchat0 <- subsetData(cellchat0)
cellchat0 <- identifyOverExpressedGenes(cellchat0)
cellchat0 <- identifyOverExpressedInteractions(cellchat0)
cellchat0 <- smoothData(cellchat0, adj = PPI.human)
cellchat0 <- computeCommunProb(cellchat0, raw.use = F)
cellchat0 <- filterCommunication(cellchat0, min.cells = 10)
cellchat0 <- computeCommunProbPathway(cellchat0)
cellchat0 <- aggregateNet(cellchat0)
cellchat0 <- netAnalysis_computeCentrality(cellchat0, slot.name = "netP")

cellchat1 <- subsetData(cellchat1)
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)
cellchat1 <- smoothData(cellchat1, adj = PPI.human)
cellchat1 <- computeCommunProb(cellchat1, raw.use = F)
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)
cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat1 <- aggregateNet(cellchat1)
cellchat1 <- netAnalysis_computeCentrality(cellchat1, slot.name = "netP")

cellchat2 <- subsetData(cellchat2)
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
cellchat2 <- smoothData(cellchat2, adj = PPI.human)
cellchat2 <- computeCommunProb(cellchat2, raw.use = F)
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)
cellchat2 <- netAnalysis_computeCentrality(cellchat2, slot.name = "netP")

save(cellchat0, file = "cellchat.processed.l2.HC.RData")
save(cellchat1, file = "cellchat.processed.l2.DCE.RData")
save(cellchat2, file = "cellchat.processed.l2.DRE.RData")

#Create merged CellChat objects to compare between diagnostic groups
object.list <- list(HC = cellchat0, DCE = cellchat1, DRE = cellchat2)
cellchatl2hc <- mergeCellChat(object.list[c(1,2)], add.names = names(object.list)[c(1,2)])
cellchatl2hr <- mergeCellChat(object.list[c(1,3)], add.names = names(object.list)[c(1,3)])
cellchatl2cr <- mergeCellChat(object.list[c(2,3)], add.names = names(object.list)[c(2,3)])
cellchatl2hcr <- mergeCellChat(object.list, add.names = names(object.list))
save(cellchatl2hc, file = "cellchat.new.l2.HC.DCE.RData")
save(cellchatl2hr, file = "cellchat.new.l2.HC.DRE.RData")
save(cellchatl2cr, file = "cellchat.new.l2.DCE.DRE.RData")
save(cellchatl2hcr, file = "cellchat.new.l2.HC.DCE.DRE.RData")