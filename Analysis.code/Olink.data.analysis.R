library(OlinkAnalyze)
library(dplyr)
library(arrow)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)

#Read in parquet file of Olink data
filepath <- "~"
olink <-  list.files(path = filepath,
                    pattern = "parquet$",
                    full.names = TRUE) |>
  lapply(OlinkAnalyze::read_NPX)  |>
  dplyr::bind_rows()

#Isolate epilepsy samples, set diagnosis and status columns
epilepsy <- olink %>% filter(PlateID == "WAEL01_Plate02_Manifest",
                            SampleQC == "PASS", AssayQC == "PASS")
clinical <- read.csv("Epilepsy.clinical.data.csv")
epilepsy$Diagnosis <- case_when(
  !epilepsy$SampleID %in% clinical$id ~ "Control",
  epilepsy$SampleID %in% ((clinical %>% filter(Status == "Medically Controlled"))$id) ~ "DCE",
  epilepsy$SampleID %in% ((clinical %>% filter(Status == "Medically Controlled"))$id) &
    !epilepsy$SampleID %in% c(25,61) ~ "DRE")

dcevshc <- epilepsy %>% filter(Diagnosis %in% c("Control", "DCE"))
drevshc <- epilepsy %>% filter(Diagnosis %in% c("Control", "DRE"))
dcevshc$Diagnosis <- factor(dcevshc$Diagnosis, levels = c("Control", "DCE"))
drevshc$Diagnosis <- factor(drevshc$Diagnosis, levels = c("Control", "DRE"))

#Run t-tests on proteins comparing epilepsy vs. control and PC vs WC
olinkdcevshc <- olink_ttest(df = dcevshc[str_detect(dcevshc$SampleID, "C")==F,], variable = 'Diagnosis')
olinkdrevshc <- olink_ttest(df = drevshc[str_detect(drevshc$SampleID, "C")==F,], variable = 'Diagnosis')

#GSEA on t-test results
dcevshc$estimate <- -1 * dcevshc$estimate
dcevshcinput = dcevshc[,"estimate"]
names(dcevshcinput) = as.character(dcevshc[,"Assay"])
dcevshcinput = sort(dcevshcinput, decreasing = TRUE)
gseadcevshc <- gseGO(geneList = dcevshcinput,
                           ont = "BP", OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",
                           exponent = 1, minGSSize = 15, maxGSSize = 200, eps = 1e-10,
                           pvalueCutoff = 0.05, pAdjustMethod = "fdr", by = "fgsea")
gseadcevshc <- gseadcevshc@result

drevshc$estimate <- -1 * drevshc$estimate
drevshcinput = drevshc[,"estimate"]
names(drevshcinput) = as.character(drevshc[,"Assay"])
drevshcinput = sort(drevshcinput, decreasing = TRUE)
gseadrevshc <- gseGO(geneList = drevshcinput,
                     ont = "BP", OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",
                     exponent = 1, minGSSize = 15, maxGSSize = 200, eps = 1e-10,
                     pvalueCutoff = 0.05, pAdjustMethod = "fdr", by = "fgsea")
gseadrevshc <- gseadrevshc@result
