library(ggplot2)
library(shiny)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(Seurat)
library(ggrepel)
library(DT)
library(ggpubr)
library(dplyr)
library(gridExtra)
library(visNetwork)
library(bslib)
library(stringr)
library(CellChat)
library(cowplot)
library(scRepertoire)
library(ggiraph)
library(reticulate)

#############
#tags$hr(style="border-color: white;")

load(file = "shiny_objects_3_1231.rdata")
load(file="MASC_plots_tables.RData")
print(dim(data_matrix))
markergenes <- rownames(data_matrix)

load("Seurat.object.with.TCR.antigen.specificity.RData")
load("Cellchat.new.Seurat.object.l2.HC.RData")
load("Cellchat.new.Seurat.object.l2.PC.RData")
load("Cellchat.new.Seurat.object.l2.WC.RData")
load("cellchat.new.l2.HC.WC.RData")
load("cellchat.new.l2.HC.PC.RData")
load("cellchat.new.l2.WC.PC.RData")
cellchatpathways <- unique(c(cellchat0@netP$pathways, cellchat1@netP$pathways,
                             cellchat2@netP$pathways, cellchatl2hc@netP$HC$pathways,
                             cellchatl2hc@netP$WC$pathways, cellchatl2hr@netP$HC$pathways,
                             cellchatl2hr@netP$PC$pathways, cellchatl2cr@netP$WC$pathways,
                             cellchatl2cr@netP$PC$pathways))
cellchatl2hc <- liftCellChat(cellchatl2hc)
cellchatl2hr <- liftCellChat(cellchatl2hr)
cellchatl2cr <- liftCellChat(cellchatl2cr, group.new = levels(cellchatl2cr@idents$WC))
cellchatcompositelist <- c(cellchatl2hc, cellchatl2hr, cellchatl2cr)
names(cellchatcompositelist) <- c("HC vs WCE", "HC vs PCE", "WCE vs PCE")
cell_meta$TCR.cloneSize <- seurat@meta.data$cloneType

#Set factor levels of input objects and variables
l1_levels <- c("B", "Mono", "CD4 T", "CD8 T", "NK", "DC", "other T", "other")
l2_levels <- c("B naive", "B intermediate", "B memory", "Plasmablast", "CD14 Mono", "CD16 Mono",
               "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "Treg", "CD8 Naive",
               "CD8 TCM", "CD8 TEM", "gdT", "MAIT", "dnT", "CD4 Proliferating",
               "CD8 Proliferating", "NK", "NK_CD56bright", "NK Proliferating",
               "cDC2", "pDC", "Platelet", "HSPC", "ILC", "ASDC")
cell_meta$predicted.celltype.l1 <- factor(cell_meta$predicted.celltype.l1,
                                          levels = l1_levels)
cell_meta$predicted.celltype.l2 <- factor(cell_meta$predicted.celltype.l2,
                                          levels = l2_levels)
seurat@meta.data$predicted.celltype.l1 <- factor(seurat@meta.data$predicted.celltype.l1,
                                                 levels = c("B", "Mono", "CD4T", "CD8T", "NK", "DC", "otherT", "other"))
seurat@meta.data$predicted.celltype.l2 <- factor(seurat@meta.data$predicted.celltype.l2,
                                                 levels = c("Bnaive", "Bintermediate", "Bmemory", "Plasmablast", "CD14Mono", "CD16Mono",
                                                            "CD4Naive", "CD4TCM", "CD4TEM", "CD4CTL", "Treg", "CD8Naive",
                                                            "CD8TCM", "CD8TEM", "gdT", "MAIT", "dnT", "CD4Proliferating",
                                                            "CD8Proliferating", "NK", "NKCD56bright", "NKProliferating",
                                                            "cDC2", "pDC", "Platelet", "HSPC", "ILC", "ASDC"))

#Edit MC to WCE and MR to PCE in metadata
cell_meta$Condition <- str_replace(cell_meta$Condition, "MC", "WCE")
cell_meta$Condition <- str_replace(cell_meta$Condition, "MR", "PCE")
l1_l2$Condition <- str_replace(l1_l2$Condition, "MC", "WCE")
l1_l2$Condition <- str_replace(l1_l2$Condition, "MR", "PCE")
lis <- list(c("HC", "WCE"), c("HC", "PCE"), c("WCE", "PCE"))
MASCinput$Group <- str_replace(MASCinput$Group, "MC", "WCE")
MASCinput$Group <- str_replace(MASCinput$Group, "MR", "PCE")
pair_df$exclude <- str_replace(pair_df$exclude, "MC", "WCE")
pair_df$exclude <- str_replace(pair_df$exclude, "MR", "PCE")
pair_df$exclude <- factor(pair_df$exclude, levels = c("HC vs WCE", "HC vs PCE", "PCE vs WCE"))
PCWCL1$Group <- str_replace(PCWCL1$Group, "MC", "WCE")
PCWCL1$Group <- str_replace(PCWCL1$Group, "MR", "PCE")
l1_l2$Condition <- str_replace(l1_l2$Condition, "MC", "WCE")
l1_l2$Condition <- str_replace(l1_l2$Condition, "MR", "PCE")
summarized_data$exclude <- str_replace(summarized_data$exclude, "MC", "WCE")
summarized_data$exclude <- str_replace(summarized_data$exclude, "MR", "PCE")
summarized_data$exclude <- factor(summarized_data$exclude,
                                  levels = c("HC vs WCE", "HC vs PCE", "PCE vs WCE"))

server <- function(input, output, session) {
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
  # vals
  
  observe({
    var <- switch(input$Variable,
                  "Cell Type (L1)" = 1,
                  "Cell Type (L2)" = 2,
                  "Sex" = 3,
                  "Condition" = 4,
                  "Seurat Clusters" = 5,
                  "TCR clonotype group" = 6)
    x <- vals[[var]]
    # Can also set the label and select items
    updateCheckboxGroupInput(session, "inCheckboxGroup",
                             label = "Variable of interest options",
                             choices = x,
                             selected = x[1]
    )
  })
  
  observe({
    var_marker_L1 <- match(input$DEG_L1, names(roc_l1_groups))
    x_marker_L1 <- roc_l1_genes_list[[var_marker_L1]]
    # Can also set the label and select items
    updateSelectInput(session, "DEG_L1_marker",
                      label = "Cell type markers",
                      choices = x_marker_L1,
                      selected = x_marker_L1[1]
    )
  })
  
  observe({
    var_marker_L2 <- match(input$DEG_L2, names(roc_l2_groups))
    x_marker_L2 <- roc_l2_genes_list[[var_marker_L2]]
    # Can also set the label and select items
    updateSelectInput(session, "DEG_L2_marker",
                      label = "Cell sub-type markers",
                      choices = x_marker_L2,
                      selected = x_marker_L2[1]
    )
  })
  
  observe({
    
    if(input$HC_vs_EP_Level =="L1"){
      x_var <- unique(hc_vs_ep_df_l1 %>% filter(p_adj <= 0.05) %>% pull(cell_type) %>% unique())
    } else{
      x_var <- unique(hc_vs_ep_df_l2 %>% filter(p_adj <= 0.05) %>% pull(cell_type) %>% unique())
    }
    # Can also set the label and select items
    updateRadioButtons(session, "HC_Epil_Cluster",
                       label = "Cell types",
                       choices = x_var,
                       selected = x_var[1]
    )
  })
  
  sub_types <- reactive({
    input$inCheckboxGroup
  })  
  sub_types_tau <- reactive({
    input$inCheckboxGroup_tau
  })
  
  sub_types_DEG_L1 <- reactive({
    input$DEG_L1_marker
  })
  sub_types_DEG_L2 <- reactive({
    input$DEG_L2_marker
  })
  
  
  ###### UMAP ###
  UMAP_plot <- function(){
    var <- switch(input$Variable,
                  "Cell Type (L1)" = "predicted.celltype.l1",
                  "Cell Type (L2)" = "predicted.celltype.l2",
                  "Sex" = "Sex",
                  "Condition" = "Condition",
                  "Seurat Clusters" = "seurat_clusters",
                  "TCR clonotype group" = "TCR.cloneSize")
    gg <- DimPlot2(group = var) + ggtitle(input$Variable)
    
    return(gg)
  }
  
  output$UMAP <- renderPlot({
    UMAP_plot()
  }, height= 500,width=600)
  
  output$downloadPlot_UMAP <- downloadHandler(
    filename = "UMAP.png",
    content = function(file) {
      png(file, height= 8, width = 9,res=1200, units = "in")
      print(UMAP_plot())
      dev.off()
    }
  )
  
  ####################  Marker genes
  
  umap_feature_L1_plot <- function(){
    plot <- FeaturePlot2(feature = sub_types_DEG_L1())
    return(plot)
  }
  
  vln_feature_L1_plot <- function(){
    #  plot <- VlnPlot(pbmc, features = sub_types_DEG_L1(), group.by="predicted.celltype.l1")+NoLegend()
    plot <- VlnPlot2(features = sub_types_DEG_L2(), group = "predicted.celltype.l1")
    
    return(plot)
  }
  
  DEG_L1_table <- function(){
    return(roc_l1 %>% filter(cluster == input$DEG_L1))
  }
  
  output$UMAP_L1 <- renderPlot({
    cell_meta$selected_column <- paste0("Not ", input$DEG_L1)
    cell_meta$selected_column[cell_meta$predicted.celltype.l1 %in% input$DEG_L1] <- input$DEG_L1
    # gg <- DimPlot(pbmc, group.by = "selected_column",label=T)+ggtitle("")+NoLegend()
    gg <- DimPlot3(group = "predicted.celltype.l1",highlight = input$DEG_L1 ,
                   position_input = "none")+ggtitle("")
    return(gg)
  })
  
  output$UMAP_feature_L1 <- renderPlot({
    umap_feature_L1_plot()
  })
  output$Vln_feature_L1 <- renderPlot({
    vln_feature_L1_plot()
  })
  output$DEG_L1_Table <- renderDataTable({
    DEG_L1_table()
  })
  
  output$downloadPlot_UMAP_Vln_feature_L1 <- downloadHandler(
    filename = paste0("UMAP_Vln_Feature",sub_types_DEG_L1(),".png"),
    content = function(file) {
      png(file, height= 8, width = 10,res=1200, units = "in")
      print(umap_feature_L1_plot() + vln_feature_L1_plot() )
      dev.off()
    }
  )
  
  output$downloadTable_ROC_results_L1 <- downloadHandler(
    filename = paste0(input$DEG_L1,"_Marker_Genes_table.csv"),
    content = function(file) {
      write.csv(DEG_L1_table(), file)
    }
  )
  #####Marker L2###
  
  umap_feature_L2_plot <- function(){
    plot <- FeaturePlot2(feature = sub_types_DEG_L2())
    return(plot)
  }
  
  vln_feature_L2_plot <- function(){
    # plot <- VlnPlot(pbmc, features = sub_types_DEG_L2(), group.by="predicted.celltype.l2")+NoLegend()
    plot <- VlnPlot2(features = sub_types_DEG_L2(), group = "predicted.celltype.l2")
    return(plot)
  }
  
  DEG_L2_table <- function(){
    return(roc_l2 %>% filter(cluster == input$DEG_L2))
  }
  
  output$UMAP_L2 <- renderPlot({
    cell_meta$selected_column <- paste0("Not ", input$DEG_L2)
    cell_meta$selected_column[cell_meta$predicted.celltype.l2 %in% input$DEG_L2] <- input$DEG_L2
    # gg <- DimPlot(pbmc, group.by = "selected_column",label=T)+ggtitle("")+NoLegend()
    gg <- DimPlot3(group = "predicted.celltype.l2", highlight = input$DEG_L2, position_input = "none")+ggtitle("")
    return(gg)
  })
  
  output$UMAP_feature_L2 <- renderPlot({
    umap_feature_L2_plot()
  })
  output$Vln_feature_L2 <- renderPlot({
    vln_feature_L2_plot()
  })
  output$DEG_L2_Table <- renderDataTable({
    DEG_L2_table()
  })
  
  output$downloadPlot_UMAP_Vln_feature_L2 <- downloadHandler(
    filename = paste0("UMAP_Vln_Feature",sub_types_DEG_L2(),".png"),
    content = function(file) {
      png(file, height= 8, width = 10,res=1200, units = "in")
      print(umap_feature_L2_plot() + vln_feature_L2_plot() )
      dev.off()
    }
  )
  
  output$downloadTable_ROC_results_L2 <- downloadHandler(
    filename = paste0(input$DEG_L2,"_Marker_Genes_table.csv"),
    content = function(file) {
      write.csv(DEG_L2_table(), file)
    }
  )
  
  ######## MASC ###
  
  #MASC_conditions_plot <- function(){
  #  if(input$Status_Comparisons == "PCE/WCE vs HC"){
  #    plot <- EpilCtrlplot
  #  }
  #  if(input$Status_Comparisons == "PCE vs WCE"){
  #    plot <- PCWCplot
  #  }
  #  return(plot)
  #}
  #MASC_conditions_table <- function(){
  #  if(input$Status_Comparisons == "PCE/WCE vs HC"){
  #    df <- MASCinput %>%
  #      mutate_if(is.numeric,function(x) round(x,4)) %>%
  #      relocate(Group)
  #  }
  #  if(input$Status_Comparisons == "PCE vs WCE"){
  #    df <- PCWCL1%>%
  #      mutate_if(is.numeric,function(x) round(x,4)) %>%
  #      dplyr::select(-Group)
  #  }
  #  return(df)
  #}
  
  #output$MASC_conditions_plot <- renderPlot({
  #  MASC_conditions_plot()
  #}, height=400,width=900)
  
  #output$MASC_conditions_results <- renderDataTable({
  #  MASC_conditions_table()
  #})
  
  #output$downloadPlot_MASC_conditions <- downloadHandler(
  #  filename = "Shiny_MASC_conditions_plot.png",
  #  content = function(file) {
  #    png(file, height= 8, width = 9,res=1200, units = "in")
  #    print(MASC_conditions_plot())
  #    dev.off()
  #  })
  #output$downloadTable_MASC_conditions <- downloadHandler(
  #  filename = "Shiny_MASC_conditions_table.csv",
  #  content = function(file) {
  #    write.csv(MASC_conditions_table(), file)
  #  })
  
  ### Treatments
  #"Valproate","Na channel blockers","Antidepressants"
  #MASC_tx_plot <- function(){
  #  if(input$Treatment == "Valproate"){
  #    plot <- MASCVPAplot
  #  }
  #  if(input$Treatment == "Na channel blockers"){
  #    plot <- MASCNAblockplot
  #  }
  #  if(input$Treatment == "Antidepressants"){
  #    plot <- MASCantideprplot
  #  }
  #  return(plot)
  #}
  #MASC_tx_table <- function(){
  #  if(input$Treatment == "Valproate"){
  #    df <- MASCVPA%>%
  #      mutate_if(is.numeric,function(x) round(x,4)) %>%
  #      relocate(Group)
  #  }
  #  if(input$Treatment == "NA channel blockers"){
  #    df <- MASCNAblock%>%
  #      mutate_if(is.numeric,function(x) round(x,4)) %>%
  #      relocate(Group)
  #  }
  #  if(input$Treatment == "Antidepressants"){
  #    df <- MASCantidepr%>%
  #      mutate_if(is.numeric,function(x) round(x,4)) %>%
  #      relocate(Group)
  #  }
  #  return(df)
  #}
  
  #output$MASC_tx_plot <- renderPlot({
  #  MASC_tx_plot()
  #}, height=400,width=900)
  
  #output$MASC_tx_results <- renderDataTable({
  #  MASC_tx_table()
  #})
  
  #output$downloadPlot_MASC_tx <- downloadHandler(
  #  filename = paste0("Shiny_MASC_tx_",input$Treatment,"_plot.png"),
  #  content = function(file) {
  #    png(file, height= 8, width = 9,res=1200, units = "in")
  #    print(MASC_conditions_plot())
  #    dev.off()
  #  })
  #output$downloadTable_MASC_tx <- downloadHandler(
  #  filename = paste0("Shiny_MASC_tx_",input$Treatment,"_table.csv"),
  #  content = function(file) {
  #    write.csv(MASC_conditions_table(), file)
  #  })
  ################
  
  sub_types_epil_ctrl <- reactive({
    input$HC_Epil_Cluster
  })
  
  volc_hc_ep_plot <- eventReactive(input$make_volcano_plot, {
    
    if(input$HC_vs_EP_Level=="L1"){
      hc_vs_ep_df <- hc_vs_ep_df_l1
    } else{
      hc_vs_ep_df <- hc_vs_ep_df_l2
    }
    colnames(hc_vs_ep_df)[grep("logfc",colnames(hc_vs_ep_df))] <- "logFC"
    plt <- hc_vs_ep_df %>%
      filter(cell_type %in% sub_types_epil_ctrl()) %>%
      mutate(differential = "no change") %>%
      mutate(differential = ifelse(p_adj <= 0.05 & logFC >0, "Up-regulated",differential)) %>%
      mutate(differential = ifelse(p_adj <= 0.05 & logFC < 0,"Down-regulated",differential)) %>%
      ggplot(aes(x=logFC, y= -log10(p_adj), color = differential))+
      ggiraph::geom_point_interactive(
        aes(tooltip = sprintf("%s\nlogFC: %s\nFDR: %s", 
                              gene, 
                              signif(logFC, digits = 3),
                              signif(p_adj, digits = 2)),
            x = logFC, y = -log10(p_adj), color = differential),
        hover_nearest = TRUE, size = 3, alpha = 0.4) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed",color ="gray")+
      theme_bw()+
      scale_color_manual(values=c("Up-regulated"="darkred", "Down-regulated"="blue","no change"="black"))+
      ggtitle(paste0(sub_types_epil_ctrl()," DEGs in epilepsy vs control"))+
      theme(legend.position = "down")
    ggiraph::girafe(ggobj = plt,
                    options = list(
                      opts_tooltip(use_fill = TRUE),
                      opts_zoom(min = 0.5, max = 5),
                      opts_sizing(rescale = TRUE),
                      opts_toolbar(saveaspng = TRUE, delay_mouseout = 2000)))
  })
  
  vln_hc_ep_plot <- function(){
    
    if(input$HC_vs_EP_Level=="L1"){
      
      genes <- hc_vs_ep_df_l1 %>%
        filter(cell_type %in% sub_types_epil_ctrl(), p_adj <= 0.05) %>%
        top_n(12, abs(logFC)) %>% pull(gene)
      #      vln_plt <- VlnPlot(subset(pbmc,subset = predicted.celltype.l1_NO_SPACES %in% sub_types_epil_ctrl()),
      #                         features = genes, ncol=as.numeric(input$Vln_Cols), group.by = "Condition")
      vln_plt <- VlnPlot2(subset_condition = cell_meta$predicted.celltype.l1_NO_SPACES %in% sub_types_epil_ctrl(),
                          features = genes, ncol_input=as.numeric(input$Vln_Cols), group = "Condition")
    } else{
      genes <- hc_vs_ep_df_l2 %>%
        filter(cell_type %in% sub_types_epil_ctrl(), p_adj <= 0.05) %>%
        top_n(15, abs(logCPM)) %>% arrange(desc(logCPM)) %>% pull(gene)
      
      #  vln_plt <- VlnPlot(subset(pbmc,subset = predicted.celltype.l2_NO_SPACES %in% sub_types_epil_ctrl()),
      #                     features = genes, ncol=as.numeric(input$Vln_Cols), group.by = "Condition")
      vln_plt <- VlnPlot2(subset_condition = cell_meta$predicted.celltype.l2_NO_SPACES %in% sub_types_epil_ctrl(),
                          features = genes, ncol_input=as.numeric(input$Vln_Cols), group = "Condition")
    }
    return(vln_plt)
  }
  output$Volcano_HC_EP <- renderGirafe({
    volc_hc_ep_plot()
  })
  output$Vln_Plots_HC_EP <- renderPlot({
    vln_hc_ep_plot()
  })
  
  output$downloadPlots_HC_vs_EP <- downloadHandler(
    filename = "Shiny_HC_vs_EP_plot.png",
    content = function(file) {
      png(file, height= 8, width = 11,res=1200, units = "in")
      print(vln_hc_ep_plot())
      dev.off()
    }
  )
  
  DEG_volc_table <- function(){
    if(input$HC_vs_EP_Level=="L1"){
      hc_vs_ep_df <- hc_vs_ep_df_l1
    } else{
      hc_vs_ep_df <- hc_vs_ep_df_l2
    }
    return(hc_vs_ep_df %>% filter(cell_type %in% sub_types_epil_ctrl()))
  }
  
  output$downloadTable_volc_results <- downloadHandler(
    filename = paste0(sub_types_epil_ctrl(),"_DGE_ep_vs_ctrl_table.csv"),
    content = function(file) {
      write.csv(DEG_volc_table(), file)
    }
  )
  
  observe({
    ft_df <- summarized_data %>% filter(exclude == input$L1_Pairwise,
                                        direction ==input$Direction_L1,
                                        level =="L1") %>% pull(cell_type)
    print(ft_df)
    updateRadioButtons(session, "DEG_L1_choice",
                       label = "Cell type",
                       choices = sort(unique(ft_df)),
                       selected =  sort(unique(ft_df))[1])
  })
  
  pair_degs_L1 <- reactive({
    input$DEG_L1_choice
  })
  
  observe({
    ft_df_l2 <- summarized_data %>% filter(exclude==input$L2_Pairwise,
                                           direction==input$Direction_L2,
                                           level =="L2") %>% pull(cell_type)
    print(ft_df_l2)
    updateRadioButtons(session, "DEG_L2_choice",
                       label = "Cell sub-type",
                       choices = sort(unique(ft_df_l2)),
                       selected =  sort(unique(ft_df_l2))[1])
  })
  
  pair_degs_L2 <- reactive({
    input$DEG_L2_choice
  })
  
  
  deg_violins_l1 <- function(){
    features_input <- features_df(input$L1_Pairwise, input$Direction_L1,"L1") %>%
      filter(cell_type == str_remove_all(pair_degs_L1()," ")) %>% arrange(desc(logCPM)) %>% pull(gene)
    # vln <- VlnPlot(subset(pbmc, predicted.celltype.l1_NO_SPACES %in% pair_degs_L1()),
    #                features = features_input[1:6], ncol=3, group.by ="Condition")
    vln <- VlnPlot2(subset_condition = cell_meta$predicted.celltype.l1_NO_SPACES %in% pair_degs_L1(),
                    features = features_input[1:6], ncol_input=3, group = "Condition")
    return(vln)
  }
  
  deg_violins_l2 <- function(){
    features_input <- features_df(input$L2_Pairwise, input$Direction_L2,"L2") %>%
      filter(cell_type == str_remove_all(pair_degs_L2()," ")) %>% arrange(desc(logCPM)) %>% pull(gene)
    # vln <- VlnPlot(subset(pbmc, predicted.celltype.l2_NO_SPACES %in% pair_degs_L2()),
    #                features = features_input[1:6], ncol=3, group.by="Condition")
    vln <- VlnPlot2(subset_condition = cell_meta$predicted.celltype.l2_NO_SPACES %in% pair_degs_L2(),
                    features = features_input[1:6], ncol_input=3, group = "Condition")
    return(vln)
  }
  deg_l1_table_pair <- function(){
    tab <- features_df(input$L1_Pairwise, input$Direction_L1,"L1")%>%
      filter(cell_type == str_remove_all(pair_degs_L1()," ")) %>%
      arrange(p_adj, desc(logCPM)) %>% dplyr::select(-any_of(c("cell_type","level","exclude","direction")))
    return(tab)
  }
  deg_l2_table_pair <- function(){
    tab <- features_df(input$L2_Pairwise, input$Direction_L2,"L2")%>%
      filter(cell_type == str_remove_all(pair_degs_L2()," ")) %>%
      arrange(p_adj, desc(logCPM)) %>% dplyr::select(-any_of(c("cell_type","level","exclude","direction")))
    return(tab)
  }
  
  output$DEG_Violins_L1 <- renderPlot({
    deg_violins_l1()
  })
  output$DEG_Violins_L2 <- renderPlot({
    deg_violins_l2()
  })
  
  output$DEG_L1_Table_Pair <- renderDataTable({
    deg_l1_table_pair() 
  })
  
  output$DEG_L2_Table_Pair <- renderDataTable({
    deg_l2_table_pair()
  })
  
  output$downloadPlot_Pairwise_L1 <- downloadHandler(
    filename = paste0("Shiny_Vln_plots",pair_degs_L1(),"_",input$L1_Pairwise,"_",input$Direction_L1,".png"),
    content = function(file) {
      png(file, height= 8, width = 11,res=1200, units = "in")
      print(deg_violins_l1())
      dev.off()
    }
  )
  output$downloadPlot_Pairwise_L2 <- downloadHandler(
    filename = paste0("Shiny_Vln_plots",pair_degs_L2(),"_",input$L2_Pairwise,"_",input$Direction_L2,".png"),
    content = function(file) {
      png(file, height= 8, width = 11,res=1200, units = "in")
      print(deg_violins_l2())
      dev.off()
    }
  )
  
  ####################
  
  dim_plot_A <- function(){
    fp <-FeaturePlot2(feature=input$marker_A)
    return(fp)
  }
  dim_plot_B <- function(){
    fp <-FeaturePlot2(feature = input$marker_B)
    return(fp)
  }
  
  vln_plot_A <- function(){
    var <- switch(input$marker_A_variable,
                  "Cell Type (L1)" = "predicted.celltype.l1",
                  "Cell Type (L2)" = "predicted.celltype.l2",
                  "Sex" = "Sex",
                  "Condition" = "Condition",
                  "Seurat Clusters" = "seurat_clusters",
                  "TCR clonotype group" = "TCR.cloneSize")
    
    vln <- VlnPlot2(features = input$marker_A, group=var)
    
    return(vln)
  }
  vln_plot_B <- function(){
    var <- switch(input$marker_B_variable,
                  "Cell Type (L1)" = "predicted.celltype.l1",
                  "Cell Type (L2)" = "predicted.celltype.l2",
                  "Sex" = "Sex",
                  "Condition" = "Condition",
                  "Seurat Clusters" = "seurat_clusters",
                  "TCR clonotype group" = "TCR.cloneSize")
    
    vln <- VlnPlot2(features = input$marker_B, group=var)
    
    return(vln)
  }
  
  output$DimPlot_A <- renderPlot({
    dim_plot_A()
  })
  output$DimPlot_B <- renderPlot({
    dim_plot_B()
  })
  
  output$Vln_Plot_A <- renderPlot({
    vln_plot_A()
  })
  output$Vln_Plot_B <- renderPlot({
    vln_plot_B()
  })
  
  
  output$dowloadPlot_Grid <- downloadHandler(
    filename = paste0(input$marker_A,"_",input$marker_B,"_grid.png"),
    content = function(file) {
      ggsave(file, (dim_plot_B()+dim_plot_B())/(vln_plot_A()+vln_plot_B()),
             height=12, width=12)
    }
  )
  
  
  #############
  UMAP_plot_2 <- function(){
    var <- switch(input$cellinfo_A_variable,
                  "Cell Type (L1)" = "predicted.celltype.l1",
                  "Cell Type (L2)" = "predicted.celltype.l2",
                  "Sex" = "Sex",
                  "Condition" = "Condition",
                  "Seurat Clusters" = "seurat_clusters",
                  "TCR clonotype group" = "TCR.cloneSize")
    gg <- DimPlot2(group = var, position_input = "bottom")
    
    return(gg)
  }
  
  UMAP_plot_3 <- function(){
    var <- switch(input$cellinfo_B_variable,
                  "Cell Type (L1)" = "predicted.celltype.l1",
                  "Cell Type (L2)" = "predicted.celltype.l2",
                  "Sex" = "Sex",
                  "Condition" = "Condition",
                  "Seurat Clusters" = "seurat_clusters",
                  "TCR clonotype group" = "TCR.cloneSize")
    gg <- DimPlot2(group = var, position_input = "bottom")
    return(gg)
  }
  
  output$Cellinfo_Plot_A <- renderPlot({
    UMAP_plot_2()
  })
  
  output$Cellinfo_Plot_B <- renderPlot({
    UMAP_plot_3()
  })
  
  
  UMAP_plot_4 <- function(){
    var <- switch(input$CxG_A_variable,
                  "Cell Type (L1)" = "predicted.celltype.l1",
                  "Cell Type (L2)" = "predicted.celltype.l2",
                  "Sex" = "Sex",
                  "Condition" = "Condition",
                  "Seurat Clusters" = "seurat_clusters",
                  "TCR clonotype group" = "TCR.cloneSize")
    gg <- DimPlot2(group = var, position_input = "bottom")
    
    return(gg)
  }
  
  output$CxG_Plot_A <- renderPlot({
    UMAP_plot_4()
  })
  
  
  dim_plot_A2 <- function(){
    fp <-FeaturePlot2(feature=input$CxG_marker_A)
    return(fp)
  }
  
  output$CxG_DimPlot_A <- renderPlot({
    dim_plot_A2()
  })
  
  ### Dot plot ###
  
  dotvar <- reactive({
    req(input$dotVar)
    switch(input$dotVar,
           "Cell Type (L1)" = "predicted.celltype.l1",
           "Cell Type (L2)" = "predicted.celltype.l2",
           "Seurat Clusters" = "seurat_clusters")
  })
  labels <- reactive({
    req(input$Genes)
    unlist(strsplit(input$Genes, ","))
  })
  dotplot <- eventReactive(input$make_plot, {
    Idents(seurat) <- dotvar()
    DotPlot(seurat, features = labels()) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())
  })
  output$DotPlot <- renderPlot({dotplot()})
  
  ### CellChat ###
  
  var <- reactive({
    req(input$composite)
  })
  pathways <- reactive({
    req(input$pathways)
  })
  
  ranknetplot <- eventReactive(input$make_ranknet_plot, {
    ranknet <- rankNet(cellchatcompositelist[[var()]],
                       mode = "comparison", measure = "weight", do.stat = TRUE,
                       signaling = pathways(), stacked = T,
                       title = paste0("Selected signaling in ", var()))
    ranknet$theme$axis.text.y$colour <- "black"
    ranknet
  })
  
  output$CellChatPlot <- renderPlot({ranknetplot()})
  
  output$CirclePlotHC <- renderPlot({
    netVisual_aggregate(cellchat0, signaling = input$pathway,
                        layout = "circle", remove.isolate = T)
  })
  
  output$CirclePlotWCE <- renderPlot({
    netVisual_aggregate(cellchat1, signaling = input$pathway,
                        layout = "circle", remove.isolate = T)
  })
  
  output$CirclePlotPCE <- renderPlot({
    netVisual_aggregate(cellchat2, signaling = input$pathway,
                        layout = "circle", remove.isolate = T)
  })
  
  comp <- reactive({
    req(input$comparison)
  })
  lrpathway <- reactive({
    req(input$lrpathway)
  })
  
  geneexpplot <- eventReactive(input$make_geneexp_plot, {
    ranknet <- rankNet(cellchatcompositelist[[comp()]],
                       mode = "comparison", measure = "weight", do.stat = TRUE,
                       signaling = lrpathway(), stacked = T,
                       title = paste0("Selected signaling in ", comp()))
    legend <- get_legend(ranknet)
    rm(ranknet)
    ccgeneexp <- plotGeneExpression(cellchatcompositelist[[comp()]],
                                    signaling = lrpathway(), split.by = "datasets",
                                    colors.ggplot = T, type = "violin")
    plot_grid(ccgeneexp,legend,rel_widths = c(9,1))
  })
  
  output$cellchatgeneexpplot <- renderPlot({geneexpplot()})
}

shinyApp(ui,server)