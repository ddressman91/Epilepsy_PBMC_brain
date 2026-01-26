library(scRepertoire)
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggsignif)
library(RColorBrewer)
library(turboGliph)
library(ggseqlogo)
library(msa)
library(Biostrings)
library(readxl)

#Load the Seurat object, contig list and the combined TCR object
setwd("~")
load("Seurat.object.with.TCR.antigen.specificity.RData")
load("Epilepsy_TCRcontigs_listby_patientID.RData")
load("Epilepsy_TCRlist_strict_TCRab.RData")

#UMAP plot with TCR clonotypes projected
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

Idents(seurat) <- "cloneType"
dimplot <- DimPlot(seurat, raster = F, cols = colorblind_vector(5)) +
  scale_color_manual(values = colorblind_vector(5), na.value="grey",
                     name = "Clone type") + theme(legend.position = "none")

#Breakdown plots of expanded clones by L2 cell type and diagnosis
l2breakdown <- clonalOccupy(seurat, x.axis = "predicted.celltype.l2")

l2plot <- l2breakdown$data %>%
  filter(cloneSize != "Single (0 < X <= 1)") %>%
  filter(n >= 10) %>%
  filter(predicted.celltype.l2 %in% c("CD4CTL", "CD4Naive", "CD4Proliferating",
                                      "CD4TCM", "CD4TEM", "CD8Naive", "CD8Proliferating",
                                      "CD8TCM", "CD8TEM", "dnT", "gdT", "MAIT", "Treg"))

l2plot <- read.csv("L2.cloneSize.breakdown.csv") #If reading in from source data

ggplot(l2plot, aes(x = predicted.celltype.l2, y = n, fill = cloneSize)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_discrete(type = c("#F0F921","#F69441","#CA4778","#7D06A5"),
                      name = "Clone type") +
  ylab("Count") + xlab("L2 cell type")  + theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1))

Idents(seurat) <- "predicted.celltype.l2"
cd8temonly <- subset(seurat, idents = "CD8TEM")
cd8tembreakdown <- clonalOccupy(cd8temonly, x.axis = "condition")

cd8templot <- cd8tembreakdown$data %>%
  filter(cloneSize != "Single (0 < X <= 1)")
cd8templot$condition <- case_when(cd8templot$condition == "0" ~ "HC",
                                  cd8templot$condition == "1" ~ "DCE",
                                  cd8templot$condition == "2" ~ "DRE")
cd8templot$condition <- factor(cd8templot$condition, levels = c("HC", "DCE", "DRE"))

cd8templot <- read.csv("CD8TEM.diagnosis.breakdown.csv") #If reading in from source data

ggplot(cd8templot, aes(x = condition, y = n, fill = cloneSize)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_discrete(type = c("#F0F921","#F69441","#CA4778","#7D06A5"),
                      name = "Clone type") +
  ylab("Count") + xlab("Condition") + theme_bw()

Idents(seurat) <- "predicted.celltype.l2"
cd4ctlonly <- subset(seurat, idents = "CD4CTL")
cd4ctlbreakdown <- clonalOccupy(cd4ctlonly, x.axis = "condition")

cd4ctlplot <- cd4ctlbreakdown$data %>%
  filter(cloneSize != "Single (0 < X <= 1)")
cd4ctlplot$condition <- case_when(cd4ctlplot$condition == "0" ~ "HC",
                                  cd4ctlplot$condition == "1" ~ "DCE",
                                  cd4ctlplot$condition == "2" ~ "DRE")
cd4ctlplot$condition <- factor(cd4ctlplot$condition, levels = c("HC", "DCE", "DRE"))

cd4ctlplot <- read.csv("CD4CTL.diagnosis.breakdown.csv") #If reading in from source data

ggplot(cd4ctlplot, aes(x = condition, y = n, fill = cloneSize)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_discrete(type = c("#F69441","#CA4778","#7D06A5"),
                      name = "Clone type") +
  ylab("Count") + xlab("Condition") + theme_bw()

#Volcano plot of large and expanded CD8+ TEM clones vs single clones
cd8temonly@meta.data$Expanded <- case_when(
  cd8temonly@meta.data$Frequency >= 21 ~ "Expanded",
  cd8temonly@meta.data$Frequency < 21 & cd8temonly@meta.data$Frequency > 1 ~ "Mid",
  cd8temonly@meta.data$Frequency == 1 ~ "Single",
  cd8temonly@meta.data$Frequency == 0 | cd8temonly@meta.data$Frequency == NA ~ "NA")

cd8tempseudo <- AggregateExpression(cd8temonly, assays = "RNA",
                                    group.by = c("PatientID", "Expanded"))
cd8tempseudo <- cd8tempseudo$RNA

cd8temmeta <- as.data.frame(table(cd8temonly@meta.data$PatientID,
                                  cd8temonly@meta.data$Expanded)) %>%
  reshape(direction = "wide", idvar = "Var1", timevar = "Var2")
colnames(cd8temmeta) <- c("PatientID", "nExpanded", "nMid", "nSingle")
cd8temmeta$Total <- rowSums(cd8temmeta[,c(2:4)])
cd8temmeta$PatientID <- as.numeric(cd8temmeta$PatientID)
cd8temmeta <- cd8temmeta %>% filter(!PatientID %in% c(25,61)) #Remove non-DRE patients
cd8temmeta <- left_join(cd8temmeta, unique(cd8temonly@meta.data %>%
                                             dplyr::select(PatientID,condition)),
                        by = "PatientID")
cd8temmeta$PatientID <- as.character(cd8temmeta$PatientID)
cd8temmeta$PatientID <- paste0("g", cd8temmeta$PatientID)
expanded <- colnames(cd8tempseudo)[str_detect(colnames(cd8tempseudo), "Expanded")==T]
expanded <- str_remove(expanded, "_.*")
datacolnames <- c(paste0(expanded, "_Expanded"), paste0(expanded, "_Single"))
deseqmeta <- cd8temmeta %>% filter(PatientID %in% expanded)
deseqdata <- cd8tempseudo[, colnames(cd8tempseudo) %in% datacolnames]
deseqmeta <- bind_rows(deseqmeta %>% mutate(Clonality = "Expanded",
                                            PatientID = paste0(PatientID, "_Expanded")),
                       deseqmeta %>% mutate(Clonality = "Single",
                                            PatientID = paste0(PatientID, "_Single")))
deseqmeta$ncells <- case_when(deseqmeta$Clonality == "Expanded" ~ deseqmeta$nExpanded,
                              deseqmeta$Clonality == "Single" ~ deseqmeta$nSingle)
deseqdata <- deseqdata[,colnames(deseqdata) %in% deseqmeta$PatientID]

deseqdata <- deseqdata[rowSums(deseqdata == 0) <= round(ncol(deseqdata)*0.8,0),]
deseqdata <- as.data.frame(deseqdata)
deseqdata$Max <- apply(X = deseqdata, MARGIN = 1, FUN = max)
deseqdata <- deseqdata %>% filter(Max >= 3)
deseqdata <- as.data.frame(deseqdata)
deseqdata <- deseqdata[,(1:(ncol(deseqdata)-1))]
deseqdata <- deseqdata[,match(deseqmeta$PatientID,colnames(deseqdata))]
deseqmeta$condition <- factor(deseqmeta$condition, levels = c(0,1,2))
deseqmeta$Clonality <- factor(deseqmeta$Clonality,
                              levels = c("Single", "Expanded"))

expvssingledds <- DESeqDataSetFromMatrix(countData = deseqdata,
                                         colData = deseqmeta,
                                         design = ~Clonality + condition + ncells)
expvssingledds <- DESeq(expvssingledds)
expvssingledds <- as.data.frame(results(expvssingledds,
                                        name = "Clonality_Expanded_vs_Single"))

expvssingledds <- read.csv("CD8TEM.expanded.vs.single.DGE.no2561.csv") #For reading in from source data

EnhancedVolcano(expvssingledds, lab = expvssingledds$gene,
                x = "log2FoldChange", y = "padj", xlim = c(-7,3.2),
                xlab = "Log2FC (expanded vs single)",
                selectLab = rownames(expvssingledds %>% filter(padj < 0.05)),
                title = "DEGs from expanded CD8+ TEM", ylim = c(0,7),
                subtitle = NULL, FCcutoff = 0.4, labSize = 3, pCutoff = 0.05,
                legendPosition = "None", caption = NULL,
                drawConnectors = T, arrowheads = F, max.overlaps = 75) +
  theme(plot.title = element_text(face = "plain"))

#Shannon entropy plots by diagnosis and seizure recency
shannon <- read_excel("Supplementary table 11.xlsx")
shannon$Diagnosios <- factor(shannon$Diagnosios, levels = c("HC", "DCE", "DRE"))

shannon$Normalized <- shannon$Shannon / max(shannon$Shannon)
shannon$mean <- case_when(
  shannon$Diagnosios == "HC" ~ mean((shannon %>% filter(Diagnosios == "HC"))$Normalized),
  shannon$Diagnosios == "DCE" ~ mean((shannon %>% filter(Diagnosios == "DCE"))$Normalized),
  shannon$Diagnosios == "DRE" ~ mean((shannon %>% filter(Diagnosios == "DRE"))$Normalized))
shannon$Diagnosios <- factor(shannon$Diagnosios, levels = c("HC", "DCE", "DRE"))

ggplot(shannon, aes(x = Diagnosis, y = Normalized)) + theme_bw() + ylim(c(0.8,1.02)) +
  geom_jitter(width = 0.2) + theme(axis.title.x = element_blank()) +
  geom_errorbar(aes(ymin = mean, ymax = mean), linewidth = 1.25, width = 0.6) +
  ylab(label = "Normalized Shannon entropy") + ggtitle("Shannon entropy") +
  geom_signif(comparisons = list(c("HC", "DRE")), annotations = "*",
              y_position = 1.01, tip_length = 0.01, vjust = 0.4, size = 1) +
  theme(axis.text.x=element_text(angle=45,hjust=1))

shannonclinical <- shannon %>% filter(!is.na(Seizure.in.past.month))

shannonclinical$Normalized <- shannonclinical$Shannon / max(shannonclinical$Shannon)
shannonclinical$mean[shannonclinical$Seizure.in.past.month == "N"] <-
  mean((shannonclinical %>% filter(Seizure.in.past.month == "N"))$Normalized)
shannonclinical$mean[shannonclinical$Seizure.in.past.month == "Y"] <-
  mean((shannonclinical %>% filter(Seizure.in.past.month == "Y"))$Normalized)
shannonclinical$Seizure.in.past.month <- factor(shannonclinical$Seizure.in.past.month,
                                           levels = c("N", "Y"))

ggplot(shannonclinical,
       aes(x = Seizure.in.past.month, y = Normalized)) + theme_bw() +
  geom_jitter(width = 0.2) + theme(axis.title.x = element_blank()) +
  geom_errorbar(aes(ymin = mean, ymax = mean), linewidth = 1.25, width = 0.8) +
  ylab(label = "Normalized Shannon entropy") + ggtitle("Shannon entropy") +
  geom_signif(comparisons = list(c("N", "Y")), annotations = "*",
              y_position = 1.005, tip_length = 0.03, vjust = 0.4, size = 1) +
  scale_x_discrete(labels = c("Last seizure > 1 mo", "Last seizure < 1 mo")) + ylim(c(0.9,1.01)) +
  theme(axis.text.x=element_text(angle=45,hjust=1))

#Plot 5 epilepsy-specific TCR clusters from GLIPH
meta <- seurat@meta.data %>% filter(!is.na(CTaa)) %>%
  dplyr::select(PatientID, condition, gender, age, predicted.celltype.l1,
                predicted.celltype.l2, CTaa, CTgene, Frequency, cloneSize)
meta$CDR3b <- sub(".*_", "", meta$CTaa)
meta$CDR3a <- sub("_.*", "", meta$CTaa)
meta$TRBV <- sub(".*_", "", meta$CTgene)
meta$TRBV <- sub("\\..*", "", meta$TRBV)
meta$TRAV <- sub("_.*", "", meta$CTgene)
meta$TRAV <- sub("\\..*", "", meta$TRAV)

tcrdb <- read.csv("C:/Users/dalli/Documents/Epilepsy data/TCR-seq analysis/TCR-db.csv", fill = T)
meta <- left_join(meta, (tcrdb %>% mutate(CDR3a = cdr3.alpha) %>%
                           dplyr::select(CDR3a, antigen.species, antigen.gene,
                                         antigen.epitope, mhc.a, mhc.b, v.alpha,
                                         j.alpha)), by = "CDR3a")
colnames(meta)[c(15:19)] <- c("CDR3a.antigen.species", "CDR3a.antigen.gene",
                              "CDR3a.antigen.epitope", "CDR3a.mhc.a", "CDR3a.mhc.b")
meta <- left_join(meta, (tcrdb %>% mutate(CDR3b = cdr3.beta) %>%
                           dplyr::select(CDR3b, antigen.species, antigen.gene,
                                         antigen.epitope, mhc.a, mhc.b, v.beta,
                                         d.beta, j.beta)), by = "CDR3b")

meta.cl <- meta %>% dplyr::select(PatientID, condition, gender, age, predicted.celltype.l1,
                                  predicted.celltype.l2, CTaa, CTgene, Frequency, cloneType,
                                  CDR3a, CDR3b, TRAV, TRBV, antigen.species, antigen.gene,
                                  antigen.epitope)
meta.cl$Paired <- "No"
meta.cl$Paired[meta.cl$CTaa %in% (metadata %>% filter(Specificity %in% c("Paired.MHCb", "Paired.MHCboth")))$CTaa] <- "Yes"
meta.cl$antigen.species[meta.cl$Paired == "No"] <- "NA"
meta.cl$antigen.gene[meta.cl$Paired == "No"] <- "NA"
meta.cl$antigen.epitope[meta.cl$Paired == "No"] <- "NA"
meta.cl <- meta.cl %>% distinct(PatientID, Frequency, CTaa, CTgene, CDR3a, CDR3b, TRAV, TRBV,
                                antigen.species, antigen.gene, antigen.epitope, predicted.celltype.l1,
                                predicted.celltype.l2, Paired)

meta.cl <- meta.cl %>% arrange(desc(Frequency))
meta.cl_avg <- meta.cl %>%
  group_by(PatientID, CDR3a, CDR3b, TRAV, TRBV) %>%
  summarise(Frequency_avg = round(mean(Frequency, na.rm = TRUE),1), .groups = 'drop') %>%
  arrange(desc(Frequency_avg))

#run turbogliph combined 
res_gliph_combined <- turboGliph::gliph_combined(cdr3_sequences = meta.cl,
                                                 min_seq_length = 0,
                                                 local_method = "fisher",
                                                 boost_local_significance = TRUE,
                                                 global_method = "fisher",
                                                 clustering_method = "GLIPH1.0",
                                                 scoring_method = "GLIPH1.0",
                                                 n_cores = 1)

selected_motifs<- res_gliph_combined$motif_enrichment$selected_motifs
selected_motifs <- selected_motifs %>%
  arrange(p.value)
selected_structs <- res_gliph_combined$global_enrichment$selected_structs
selected_structs <- selected_structs %>%
  arrange(p.value)
connections <- res_gliph_combined$connections
cluster_properties <- res_gliph_combined$cluster_properties
cluster_properties <- cluster_properties %>%
  arrange(total.score)
cluster_list <- res_gliph_combined$cluster_list
ordered_cluster_list <- cluster_list[match(cluster_properties$tag, names(cluster_list))]

generate_motif <- function(df) {
  sequences <- df$CDR3b  # Extract the CDR3b column from the dataframe
  patients <- unique(df$patient)
  
  if (length(sequences) == 1) {
    return(list(motif = sequences, patients = paste(patients, collapse = ", ")))  # If there's only one sequence, return it as the motif
  }
  
  # Convert the sequences to AAStringSet (assuming amino acid sequences)
  sequences_set <- AAStringSet(sequences)
  
  # Perform multiple sequence alignment (MSA)
  alignment <- msa(sequences_set)
  
  # Generate the consensus sequence (motif)
  consensus_seq <- consensusString(alignment, threshold = 1)  # Use a 100% agreement threshold
  
  consensus_seq <- gsub("\\?", "X", consensus_seq)
  
  
  return(list(motif = consensus_seq, patients = paste(patients, collapse = ", ")))
}

motif_data <- lapply(ordered_cluster_list, generate_motif)

motifs_df <- data.frame(
  Cluster = names(motif_data),
  Motif = sapply(motif_data, function(x) x$motif),
  Patients = sapply(motif_data, function(x) x$patients),
  stringsAsFactors = FALSE
)

row.names(motifs_df) <- NULL
motifs_df$Cluster <- gsub("CRG-","",motifs_df$Cluster)

names(ordered_cluster_list) <- gsub("CRG-","",names(ordered_cluster_list))
motifs_df <- motifs_df[match(names(ordered_cluster_list), motifs_df$Cluster), ]

# Function to isolate antigen-specific TCRs in the cluster list and create sequence logos
specifictcrs <- unique((meta.cl %>% filter(Paired == "Yes"))$CTaa)
specificbetas <- str_remove(specifictcrs, ".*_")
specificlist <- ordered_cluster_list[names(ordered_cluster_list) %in% specificbetas]
specificlist <- specificlist[c(1,4,5,8,10)]
specificmotifs <- motifs_df[match(names(specificlist), motifs_df$Cluster), ]
plot_sequence_logos <- function(specificlist) {
  
  cluster_names <- names(specificlist)
  cluster_motif_map <- setNames(specificmotifs$Motif, cluster_names)
  head(cluster_motif_map)
  
  # Prepare a list to store the sequences for each cluster
  seq_data <- list()
  cluster_motif_map <- setNames(specificmotifs$Motif, cluster_names)
  
  # Loop through the first 'n_clusters' and extract the 'CDR3b' sequences
  for (i in 1:length(specificlist)) {
    cluster <- specificlist[[i]]  # Get the ith cluster data frame
    seq_data[[paste("Cluster", i,": ",specificmotifs$Motif[i])]] <- cluster$CDR3b  # Store the CDR3b sequences, named by cluster
  }
  
  # Generate the sequence logos for the clusters
  ggseqlogo(seq_data, method = "bits") +
    facet_wrap(~seq_group, scales = "free", ncol = 1) +
    theme_minimal() + 
    ggtitle("CDR3b Motifs of shared epilepsy-specific TCR clusters") 
}

plot_sequence_logos(specificlist)
