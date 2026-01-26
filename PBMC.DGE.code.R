library(tidyverse)
library(edgeR)
library(Seurat)
library(SummarizedExperiment)
library(patchwork)
library(SingleCellExperiment)
library(Matrix)
library(Matrix.utils)
library(limma)

#### Pseudobulking function
convert_data_to_pseudobulk <- function(seurat_object, sample_variable, #pseudo-individual level
                                        cluster_variable, # clusters
                                        subset_variable=NULL, subset_include=NULL ){
  if(!is.null(subset_variable) && !is.null(subset_include)){
  eval(parse(text = paste0("seurat <- subset(seurat_object, subset = ",subset_variable," %in% c('", paste0(subset_include,collapse="','"),"'))" )))
  } else{
    seurat <- seurat_object
  }
  eval(parse(text= paste0("seurat <- subset(seurat, subset = ", cluster_variable," %in% c(NA), invert=T)")))

keep <- table(seurat@meta.data[[cluster_variable]])[table(seurat@meta.data[[cluster_variable]]) > 75] %>% names()
seurat <- seurat[,seurat@meta.data[[cluster_variable]] %in% keep] 
rownames(seurat@meta.data) <- colnames(seurat) 
print(keep)
  Idents(object = seurat) <- cluster_variable
  
  counts <- GetAssayData(object = seurat, layer = "counts", assay="SCT")
  metadata <- seurat@meta.data

  if(is.factor(Idents(seurat))){
    metadata$cluster_id <- seurat@active.ident
    
  } else{
  metadata$cluster_id <- factor(seurat@active.ident)
  }
  # Create single cell experiment object
  sce <- SingleCellExperiment(assays = list(counts = counts),colData = metadata)
  sce <- sce[rowSums(counts(sce) > 0) >= 30, ]  #### keep genes present in at least 30 cells
  kids <- purrr::set_names(levels(sce$cluster_id))
  nk <- length(kids)
  colData(sce)[[sample_variable]] <- paste0(colData(sce)[[sample_variable]],"_",sce$cluster_id)
  # Named vector of sample names
  sids <- purrr::set_names(levels(as.factor(colData(sce)[[sample_variable]])))
  sce$sampleInfo <- factor(colData(sce)[[sample_variable]],
                           labels = unique(colData(sce)[[sample_variable]]),
                           levels= unique(colData(sce)[[sample_variable]]))
  ns <- length(sids)

  n_cells <- table(sce$sampleInfo) %>%  as.vector()
  names(n_cells) <- names(table(sce$sampleInfo))
  
  ## Match the named vector with metadata to combine it
  m <- match(names(n_cells), sce$sampleInfo)
  ei <- data.frame(colData(sce)[m, ], 
                   n_cells, row.names = NULL) %>% 
    dplyr::select("sampleInfo", "cluster_id", "n_cells")
  # Aggregate the counts per sample_id and cluster_id
  
  # Subset metadata to only include the cluster and sample IDs to aggregate across
  groups <- colData(sce)[, c("cluster_id", "sampleInfo")]
  groups$sampleInfo <- factor(groups$sampleInfo)

  pb <- aggregate.Matrix(t(counts(sce)), 
                         groupings = groups, fun = "sum") 
  
  # Not every cluster is present in all samples; create a vector that represents how to split samples
  splitf <- sapply(stringr::str_split(rownames(pb), 
                                      pattern = "_",n = 2), `[`, 1)
 pb <- split.data.frame(pb,factor(splitf)) %>%
    lapply(function(u) 
      magrittr::set_colnames(t(u), str_extract(rownames(u),paste0(unique(ei$sampleInfo),collapse="|"))))
pb <- lapply(pb, function(u)
   u[rowSums(u > 1) >= 10, ] #### Keep genes with more than a single count across ten pseudo samples after pseudobulking
) 
  output <- list(pb = pb, ei = ei)
  return(output)
}

#### Combine pseudobulk and metadata into single object 
add_metadata <- function(d2p_output, metadata, metadata_sample_var,
                         var_of_interest, ctrl_vars=NULL, var_of_interest_reference){
  pb <- d2p_output[["pb"]]
  ei <- d2p_output[["ei"]]
  ei <- dp2$ei %>%
    as_tibble() %>%
    mutate(sampleInfo = factor(sampleInfo)) %>%
    left_join(metadata %>% dplyr::select(any_of(c(metadata_sample_var,var_of_interest,ctrl_vars))) %>% distinct(),
 by = c("sampleInfo" = metadata_sample_var)) %>%
    dplyr::select(any_of(c("sampleInfo",var_of_interest, ctrl_vars))) %>%
    distinct()
  ei[[var_of_interest]] <- factor(ei[[var_of_interest]], levels = c(var_of_interest_reference, unique(ei[[var_of_interest]])[!(unique(ei[[var_of_interest]]) %in% var_of_interest_reference)]),
                          labels= c(0,1)) %>% as.character() %>% as.numeric()
  output <- list(pb = pb,ei = ei)
  
}

##### Performs the limma- function
limma_fun <- function(pb, cluster_class_specific, var_of_interest,ctrl_var){

colnames(pb$pb[[cluster_class_specific]]) <- str_remove(colnames(pb$pb[[cluster_class_specific]]),"^JM")
d0 <- DGEList(pb$pb[[cluster_class_specific]])

print("NormFact")
d0 <- calcNormFactors(d0)
cutoff <- 1
d <- d0
meta <- pb$ei %>% mutate(sampleInfo = paste0(as.character(sampleInfo),"_",cluster_class_specific)) %>%
 filter(sampleInfo %in% colnames(pb$pb[[cluster_class_specific]]))
meta <- meta[order(match(meta$sampleInfo,colnames(pb$pb[[cluster_class_specific]]))),]
form <- as.formula(paste0("~1+",paste0("meta$",c(var_of_interest,ctrl_var),collapse="+")))
if(length(unique(meta[[var_of_interest]]))>1){
mm <- model.matrix(form)

y <- voom(d, mm, plot = F)
print("lmfit")
fit <- lmFit(y, mm)
tmp <- eBayes(fit)
top.table <- topTable(tmp, sort.by = "none", n = Inf,coef = paste0("meta$",var_of_interest)) %>%
 rownames_to_column("gene") %>% mutate(cluster=cluster_class_specific)
return(top.table)
}
}

#####################

#Data loading and brief pre-processing to handle potential character errors

directory <- "~"
load(file = "Seurat.object.with.TCR.antigen.specificity.RData")

pbmc <- seurat
rm(seurat)
pbmc$predicted.celltype.l1 <- str_remove_all(pbmc$predicted.celltype.l1, " |_")
pbmc$predicted.celltype.l2 <- str_remove_all(pbmc$predicted.celltype.l2, " |_")

#Epil_Ctrl column parsing for comparing epilepsy to control
pbmc@meta.data$Epil_Ctrl <- NA
pbmc@meta.data$Epil_Ctrl[pbmc@meta.data$condition == 1] <- 1 #Run for DCE vs HC
pbmc@meta.data$Epil_Ctrl[pbmc@meta.data$condition == 2 &
                           !(pbmc@meta.data$PatientID %in% c(25,61))] <- 1 #Run for DRE vs HC and DRE vs DCE
pbmc@meta.data$Epil_Ctrl[pbmc@meta.data$condition == 1] <- 0 #Run for DRE vs DCE
pbmc@meta.data$Epil_Ctrl[pbmc@meta.data$condition == 0] <- 0 #Run for DRE or DCE vs HC

#Epil_Ctrl column parsing for interaction term between sex and diagnosis
pbmc@meta.data$Diagnosis <- NA
pbmc@meta.data$Epil_Ctrl[pbmc@meta.data$condition == 1] <- 1 #Run for DCE vs HC
pbmc@meta.data$Diagnosis[pbmc@meta.data$condition == 2 &
                           !(pbmc@meta.data$PatientID %in% c(25,61))] <- 1 #Run for DRE vs HC and DRE vs DCE
pbmc@meta.data$Diagnosis[pbmc@meta.data$condition == 1] <- 0 #Run for DRE vs DCE
pbmc@meta.data$Diagnosis[pbmc@meta.data$condition == 0] <- 0 #Run for DRE or DCE vs HC
pbmc@meta.data$Epil_Ctrl <- pbmc@meta.data$Diagnosis * pbmc@meta.data$Sex

pbmc$predicted.celltype.l2[pbmc$predicted.celltype.l2=="NK_CD56bright"] <- "NKCD56bright"
pbmc$predicted.celltype.l1[pbmc$predicted.celltype.l1=="other"] <- "undetermined"
pbmc$predicted.celltype.l2[pbmc$predicted.celltype.l2=="NKProliferating"] <- "Proliferating"
pbmc$predicted.celltype.l2[pbmc$predicted.celltype.l2=="NK"] <- "NKCD56dim"

######################

cats <- c("l1","l2")

##### DCE vs HC
for(cat in cats){
  sub <- subset(pbmc, subset = !is.na(pbmc@meta.data$Epil_Ctrl))
  dp2 <- convert_data_to_pseudobulk(
    sub, sample_variable = "PatientID",
    cluster_variable =paste0("predicted.celltype.", cat))
  meta <- sub@meta.data %>% mutate(PatientID = factor(PatientID))
  print("Metadata")
  
  dp2$ei$sampleInfo <- parse_number(str_extract(dp2$ei$sampleInfo,"^.+_"))
  dp2 <- add_metadata(d2p_output = dp2, metadata = meta,
                    metadata_sample_var = "PatientID",
                    var_of_interest="Epil_Ctrl", var_of_interest_reference= 0,
                    ctrl_vars = c("Age","Sex"))
  print("voom")

results <- map_dfr(names(dp2$pb), function(x) {
print(x)
if(ncol(dp2$pb[[x]])>=5){
limma_fun(dp2, cluster_class_specific=x,var_of_interest = "Epil_Ctrl",ctrl_var =c("Age","Sex"))
}
}
)
results %>% write.csv(paste0(directory,"limma_voom_",cat,"_condition1_vs_condition0.csv"), row.names = F)

}