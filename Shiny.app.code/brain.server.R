server <- function(input, output, session) {
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
  # vals
  
  observe({
    var <- switch(input$BrainVar,
                   "Cell Type" = "CellType",
                   "Region" = "Region",
                   "Diagnosis" = "Diagnosis",
                   "TCR clonotype group" = "cloneSize")
  })
  
  #reactives
  
  sub_types_tle_ctrl <- reactive({
    input$HC_TLE_Cluster
  })
  
  ###### Brain UMAP
  UMAP_brain <- function(){
    var <- switch(input$BrainVar,
                  "Cell Type" = "CellType",
                  "Region" = "Region",
                  "Diagnosis" = "Diagnosis",
                  "TCR clonotype group" = "cloneSize")
    if(input$BrainVar=="cloneSize"){
      gg <- DimPlot2(group = var) +
        scale_color_manual(values = colorblind_vector(5), na.value="grey", name = "Clone type")
        ggtitle(input$BrainVar)
    } else {
      gg <- DimPlot2(group = var) + ggtitle(input$BrainVar)
    }
    
    return(gg)
  }
  
  output$UMAP <- renderPlot({
    UMAP_brain()
  }, height= 500,width=800)
  
  output$downloadPlot_UMAP <- downloadHandler(
    filename = "UMAP.png",
    content = function(file) {
      png(file, height= 8, width = 9,res=1200, units = "in")
      print(UMAP_brain())
      dev.off()
    }
  )
  
  ##### Brain cell x gene plot
  
  UMAP_brain_2 <- function(){
    var <- switch(input$CxG_Brain_UMAP,
                  "Cell Type" = "CellType",
                  "Region" = "Region",
                  "Diagnosis" = "Diagnosis",
                  "TCR clonotype group" = "cloneSize")
    if(input$BrainVar=="cloneSize"){
      gg <- DimPlot2(group = var)
      scale_color_manual(values = colorblind_vector(5), na.value="grey", name = "Clone type") +
        ggtitle(input$CxG_Brain_UMAP)
    } else {
      gg <- DimPlot2(group = var) + ggtitle(input$CxG_Brain_UMAP)
    }
    
    return(gg)
  }
  
  output$CxG_Brain_A <- renderPlot({
    UMAP_brain_2()
  })
  
  
  Feature_plot_Brain <- function(){
    fp <-FeaturePlot2(feature = input$CxG_marker_Brain)
    return(fp)
  }
  
  output$CxG_FeaturePlot_Brain <- renderPlot({
    Feature_plot_Brain()
  })
  
  ####################  Brain marker genes
  
  Vln_Brain <- function(){
    plot <- VlnPlot2(features = input$Vln_Brain_marker, group = "CellType") +
      NoLegend()
    return(plot)
  }
  
  output$Vln_feature_Brain <- renderPlot({
    Vln_Brain()
  })
  
  DEG_Brain_cluster_table <- function(){
    return(mgclust %>% filter(cluster == input$Vln_MG_celltype, p_val_adj < 0.05) %>%
             arrange(desc(avg_log2FC)))
  }
  
  output$DEG_MG_Table <- renderDataTable({
    DEG_Brain_cluster_table()
  })
  
  output$downloadTable_VlnPlot_Brain <- downloadHandler(
    filename = paste0(input$DEG_L1,"_Marker_Genes_table.csv"),
    content = function(file) {
      write.csv(DEG_Brain_cluster_table(), file)
    }
  )
  
  ####### Brain DGE
  
  volc_hc_tle_plot <- eventReactive(input$make_volcano_plot, {
    
    plt <- braincellsdds[[sub_types_tle_ctrl()]] %>%
      mutate(differential = "no change") %>%
      mutate(differential = ifelse(padj <= 0.05 & log2FoldChange >0, "Up-regulated",differential)) %>%
      mutate(differential = ifelse(padj <= 0.05 & log2FoldChange < 0,"Down-regulated",differential)) %>%
      ggplot(aes(x=log2FoldChange, y= -log10(padj), color = differential))+
      ggiraph::geom_point_interactive(
        aes(tooltip = sprintf("%s\nlogFC: %s\nFDR: %s", 
                              gene, 
                              signif(log2FoldChange, digits = 3),
                              signif(padj, digits = 2)),
            x = log2FoldChange, y = -log10(padj), color = differential),
        hover_nearest = TRUE, size = 3, alpha = 0.4) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed",color ="gray")+
      theme_bw()+
      scale_color_manual(values=c("Up-regulated"="darkred", "Down-regulated"="blue","no change"="black"))+
      ggtitle(paste0("Differential genes in DRE vs control brain for ", sub_types_tle_ctrl()))+
      theme(legend.position = "down")
    ggiraph::girafe(ggobj = plt,
                    options = list(
                      opts_tooltip(use_fill = TRUE),
                      opts_zoom(min = 0.5, max = 5),
                      opts_sizing(rescale = TRUE),
                      opts_toolbar(saveaspng = TRUE, delay_mouseout = 2000)))
  })
  
  vln_hc_tle_plot <- function(){
    
    genes <- braincellsdds[[sub_types_tle_ctrl()]] %>%
      filter(padj < 0.05) %>% top_n(15, abs(log2FoldChange)) %>%
      arrange(desc(log2FoldChange)) %>% pull(gene)
      
      vln_plt <- VlnPlot2(subset_condition = cell_meta$CellType %in% sub_types_tle_ctrl(),
                         features = genes, ncol_input = 5, group = "Diagnosis")
      return(vln_plt)
      }

  output$Volcano_HC_TLE <- renderGirafe({
    volc_hc_tle_plot()
  })
  output$Vln_Plots_HC_TLE <- renderPlot({
    vln_hc_tle_plot()
  })
  
  output$downloadPlots_HC_vs_TLE <- downloadHandler(
    filename = paste0("Shiny_HC_vs_DRE_", input$HC_TLE_Cluster,
                      "_violin_plot.png"),
    content = function(file) {
      png(file, height= 8, width = 11,res=1200, units = "in")
      print(vln_hc_tle_plot())
      dev.off()
    }
  )
  
  DEG_volc_table <- function(){
    return(braincellsdds[[sub_types_tle_ctrl()]])
  }
  
  output$downloadTable_volc_results <- downloadHandler(
    filename = paste0(sub_types_tle_ctrl(),"_DGE_DRE_vs_control_table.csv"),
    content = function(file) {
      write.csv(DEG_volc_table(), file)
    }
  )
  
  ### CellChat ###
  
  pathways <- reactive({
    req(input$pathways)
  })
  
  ranknetplot <- eventReactive(input$make_ranknet_plot, {
    ranknet <- rankNet(cellchatepctrl,
                       mode = "comparison", measure = "weight", do.stat = TRUE,
                       signaling = pathways(), stacked = T,
                       title = "Selected signaling in DRE vs Control")
    ranknet$theme$axis.text.y$colour <- "black"
    ranknet
  })
  
  output$BrainCellChatPlot <- renderPlot({ranknetplot()})
  
  output$BrainCirclePlotCtrl <- renderPlot({
    netVisual_aggregate(cellchatctrl, signaling = input$pathway,
                        layout = "circle", remove.isolate = T)
  })
  
  output$BrainCirclePlotTLE <- renderPlot({
    netVisual_aggregate(cellchatepilepsy, signaling = input$pathway,
                        layout = "circle", remove.isolate = T)
  })
  
  lrpathway <- reactive({
    req(input$brainlrpathway)
  })
  
  geneexpplot <- eventReactive(input$make_geneexp_plot, {
    ranknet <- rankNet(cellchatepctrl,
                       mode = "comparison", measure = "weight", do.stat = TRUE,
                       signaling = lrpathway(), stacked = T,
                       title = "Selected signaling in DRE vs Control")
    legend <- get_legend(ranknet)
    rm(ranknet)
    ccgeneexp <- plotGeneExpression(cellchatepctrl,
                                    signaling = lrpathway(), split.by = "datasets",
                                    colors.ggplot = T, type = "violin")
    plot_grid(ccgeneexp,legend,rel_widths = c(9,1))
  })
  
  output$braincellchatgeneexpplot <- renderPlot({geneexpplot()})
}
