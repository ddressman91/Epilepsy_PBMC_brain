library(Seurat)
library(scRepertoire)
library(dplyr)
library(stringr)

#Read in filtered_contig_associations files from each batch
batches <- c("WHE001", "WHE002", "WHE004", "WHE005", "WHE006", "WHE007", 
             "WHE008", "WHE009", "WHE010", "WHE011", "WHE093", "WHEC001", 
             "WHEC002", "WHEC003")
contigs <- vector("list", 14)
names(contigs) <- batches
tcrdf <- c()
for (batch in batches) {
  contigs[[batch]] <- read.csv(paste0("~/", batch, 
                                      "_filtered_contig_annotations.csv", fill = T)) %>%
    mutate(barcode = paste0(batch, "_", barcode), Batch = batch)
  tcrdf <- bind_rows(tcrdf, contigs[[batch]])
}

#Create TCR contig list, add diagnosis data, combine contigs
load("EpilepsyPBMC_RObject.rda")
contiglist <- createHTOContigList(contig = tcrdf, sc = pbmc,
                                  group.by = "PatientID")

conditions <- pbmc@meta.data[!duplicated(pbmc@meta.data$PatientID), "condition"]
conditions <- str_replace(conditions, "0", "HC")
conditions <- str_replace(conditions, "1", "DCE")
conditions <- str_replace(conditions, "2", "DRE")

combined <- combineTCR(contiglist, samples = names(contiglist), ID = conditions,
                       removeNA = T, removeMulti = T)

save(contiglist, file = "Epilepsy_TCRcontigs_listby_patientID.RData")
save(combined, file = "Epilepsy_TCRlist_strict_TCRab.RData")

#Calculate Shannon entropy index to measure TCR diversity
combinedtrim <- combined[c(1:33,35:62,64:85,87:91,93,95:97,99:106)] #Remove samples with low T cell counts
diversityplot <- clonalDiversity(combinedtrim, cloneCall = "strict",
                               group.by = "sample", x.axis = "ID", n.boots = 100)
shannon <- diversityplot$data