ui <- fluidPage(
  theme = shinythemes::shinytheme("united"),
  titlePanel("PBMC data in epileptic patients and controls"),
  
  tabsetPanel(               
    tabPanel("PBMC Data Overview", h3("Select cluster"),
             sidebarLayout(
               sidebarPanel(
                 radioButtons("Variable", h4("Select variable to visualize"),
                              choices = list("Cell Type (L1)","Cell Type (L2)",
                                             "Sex", "Condition",
                                             "Seurat Clusters","TCR clonotype group"),
                              selected = "Cell Type (L1)"),
                 #checkboxGroupInput("inCheckboxGroup", "Input checkbox",
                 #                    unlist(vals), selected = unlist(vals)[1]),
                 downloadButton('downloadPlot_UMAP', 'Download Plot',size="s")
                 
               ),
               mainPanel(plotOutput("UMAP")),
             )
             
    ),
    tabPanel("Explore markers",
             fluidRow(
               column(6,
                      selectizeInput("marker_A",h3("Gene A"),choices = markergenes,selected = "CD8A"),
                      # textInput("marker_A",h3("Gene A"),value = "CD8A"),
                      
                      plotOutput("DimPlot_A"),
                      hr(),
                      selectInput("marker_A_variable",h4("Variable of interest"),
                                  choices= c("Cell Type (L1)","Cell Type (L2)",
                                             "Sex", "Condition",
                                             "Seurat Clusters","TCR clonotype group"),
                                  selected = "Condition"),
                      plotOutput("Vln_Plot_A"),
                      hr(),
                      downloadButton("dowloadPlot_Grid","Download Figures",size="s")),
               column(6,
                      selectizeInput("marker_B",h3("Gene B"),choices = markergenes,selected = "CD79A"),
                      #  textInput("marker_B",h3("Gene B"),value = "CD79A"),
                      plotOutput("DimPlot_B"),
                      hr(),
                      selectInput("marker_B_variable",h4("Variable of interest"),
                                  choices =c("Cell Type (L1)","Cell Type (L2)",
                                             "Sex", "Condition",
                                             "Seurat Clusters","TCR clonotype group"),
                                  selected = "Condition" ),
                      plotOutput("Vln_Plot_B"))
             )),
    tabPanel("Explore cell information",
             fluidRow(
               column(6,
                      selectInput("cellinfo_A_variable",h4("Variable of interest"),
                                  choices= c("Cell Type (L1)","Cell Type (L2)",
                                             "Sex", "Condition",
                                             "Seurat Clusters","TCR clonotype group"),
                                  selected = "Condition"),
                      plotOutput("Cellinfo_Plot_A"),
                      hr()),
               column(6,
                      selectInput("cellinfo_B_variable",h4("Variable of interest"),
                                  choices =c("Cell Type (L1)","Cell Type (L2)",
                                             "Sex", "Condition",
                                             "Seurat Clusters","TCR clonotype group"),
                                  selected = "Cell Type (L1)" ),
                      plotOutput("Cellinfo_Plot_B")),
               hr()
             )),
    tabPanel("Explore cell x gene information",
             fluidRow(
               column(6,
                      selectInput("CxG_A_variable",h3("Variable of interest"),
                                  choices= c("Cell Type (L1)","Cell Type (L2)",
                                             "Sex", "Condition",
                                             "Seurat Clusters","TCR clonotype group"),
                                  selected = "Condition"),
                      plotOutput("CxG_Plot_A")),
               column(6,
                      selectizeInput("CxG_marker_A",h3("Gene input"),choices = markergenes,selected = "CD8A"),
                      #textInput("CxG_marker_A",h3("Gene input"),value = "CD8A"),
                      plotOutput("CxG_DimPlot_A")))),
    tabPanel("Marker genes",
             h3("Cell Type (Level 1)"),
             fluidRow(
               column(2,
                      selectInput("DEG_L1", h4("Cell type"),
                                  choices = lapply(names(roc_l1_groups),function(x) x), selected = 1),
                      hr(),
                      selectInput("DEG_L1_marker", h4("Cell type marker"),
                                  choices = lapply(names(roc_l1_genes),function(x) x)),
                      hr(),
                      downloadButton('downloadPlot_UMAP_Vln_feature_L1', 'Download Plots',size="s"),
                      hr(),
                      downloadButton('downloadTable_ROC_results_L1', 'Download Table',size="s")
                      
               ),
               column(5,plotOutput("UMAP_L1")
               ),
               column(5,plotOutput("UMAP_feature_L1")
               )),
             hr(),
             fluidRow(column(12, plotOutput("Vln_feature_L1"))),
             hr(),
             dataTableOutput("DEG_L1_Table"),
             hr(),
             h3("Cell Type (Level 2)"),
             fluidRow(
               column(2,
                      selectInput("DEG_L2", h4("Cell sub-type"),
                                  choices = lapply(names(roc_l2_groups),function(x) x), selected = 1),
                      hr(),
                      selectInput("DEG_L2_marker", h4("Cell sub-type marker"),
                                  choices = lapply(names(roc_l2_genes),function(x) x)),
                      hr(),
                      downloadButton('downloadPlot_UMAP_Vln_feature_L2', 'Download Plots',size="s"),
                      hr(),
                      downloadButton('downloadTable_ROC_results_L2', 'Download Table',size="s")
                      
               ),
               column(5,plotOutput("UMAP_L2")
               ),
               column(5,plotOutput("UMAP_feature_L2")
               )),hr(),
             fluidRow(column(12, plotOutput("Vln_feature_L2"))),
             hr(),
             dataTableOutput("DEG_L2_Table")
    ),
    #tabPanel("Dot plot of selected markers",
    #         sidebarLayout(
    #           sidebarPanel(
    #             radioButtons("dotVar", h4("Select variable to visualize"),
    #                          choices = list("Cell Type (L1)","Cell Type (L2)",
    #                                         "Seurat Clusters"),
    #                          selected = "Cell Type (L1)"),
    #             textInput("Genes", h4("Type list of genes to display, separated by commas"),
    #                       value = "CD3E"),
    #             actionButton("make_plot", label = "Create dot plot")),
    #           mainPanel(plotOutput("DotPlot"))
    #         )),
    #tabPanel("MASC Analysis",
    #         h3("Cell type enrichment between healthy controls, medically controlled, and poorly controlled patients"),
    #         fluidRow(column(12,plotOutput("MASC_conditions_plot"))),
    #         hr(),
    #         fluidRow(
    #           column(3,
    #                  radioButtons("Status_Comparisons",h4("Status comparisons"),
    #                               choices = list("PCE/WCE vs HC","PCE vs WCE"),
    #                               selected = "PCE/WCE vs HC"),hr(),
    #                  downloadButton('downloadPlot_MASC_conditions', 'Download Plot',size="s"),hr(),
    #                  downloadButton('downloadTable_MASC_conditions', 'Download Table',size="s")
    #           ),
    #           column(9,dataTableOutput("MASC_conditions_results"))),
    #         hr(),
    #         h3("Cell type enrichment across patients by medication selection"),
    #         fluidRow(column(12,plotOutput("MASC_tx_plot"))),
    #         hr(),
    #         fluidRow(
    #           column(3,
    #                  radioButtons("Treatment",h4("Medication selection"),
    #                               choices = list("Valproate","Na channel blockers","Antidepressants"),
    #                               selected = "Valproate"),hr(),
    #                  downloadButton('downloadPlot_MASC_tx', 'Download Plot',size="s"),hr(),
    #                  downloadButton('downloadTable_MASC_tx', 'Download Table',size="s")
    #           ),
    #           column(9,dataTableOutput("MASC_tx_results")))),
    tabPanel("Gene expression changes in epilepsy",
             h3("Select major cell type"),
             fluidRow(
               column(2,
                      radioButtons("HC_vs_EP_Level",label = h4("Cell Type Level"),
                                   choices = c("L1","L2"), selected = "L1"), hr(),
                      radioButtons("HC_Epil_Cluster",label="Cell type",
                                   choices=lapply(names(roc_l1_groups),function(x) x) , selected = 1),
                      hr(),
                      actionButton("make_volcano_plot", label = "Create volcano plot"),
                      downloadButton('downloadPlots_HC_vs_EP', 'Download violin plot',size="s"),
                      downloadButton('downloadTable_volc_results', 'Download Table',size="s")
               ),
               column(10,girafeOutput("Volcano_HC_EP")
               )),
             fluidRow(
               column(2,
                      radioButtons("Vln_Cols", h4("Number of columns"),
                                   choices= lapply(1:5, function(x) x), selected=4)),
               column(10,plotOutput("Vln_Plots_HC_EP"))
               
             )),
    tabPanel("Group pairwise comparisons",
             h3("Broad cell type"),
             fluidRow(
               column(2,
                      radioButtons("L1_Pairwise", label = h4("Comparison"),
                                   choices = c("HC vs DCE","HC vs DRE","DRE vs DCE"),
                                   selected = "HC vs DRE"),
                      hr(),
                      radioButtons("Direction_L1",h4("Direction"),
                                   choices = c("Up-regulated","Down-regulated"), 
                                   selected = "Down-regulated"),
                      hr(),
                      radioButtons("DEG_L1_choice", h4("Cell Type"),
                                   choices= l1_names, 
                                   selected = "B"),
                      hr(),
                      downloadButton("downloadPlot_Pairwise_L1",'Download Plot',size="s")),
               column(10,
                      plotOutput("DEG_Violins_L1"),
                      dataTableOutput("DEG_L1_Table_Pair"))
             ),
             hr(),
             h3("Sub cell type"),
             fluidRow(
               column(2, 
                      radioButtons("L2_Pairwise", label = h4("Comparison"),
                                   choices = c("HC vs DCE","HC vs DRE","DRE vs DCE"),
                                   selected = "HC vs DRE"),
                      hr(),
                      radioButtons("Direction_L2",h4("Direction"),
                                   choices = c("Up-regulated","Down-regulated"), 
                                   selected = "Down-regulated"),
                      hr(),
                      radioButtons("DEG_L2_choice", h4("Sub Cell Type"),
                                   choices= l2_names, 
                                   selected = "Bintermediate"),
                      hr(),
                      downloadButton("downloadPlot_Pairwise_L2",'Download Plot',size="s")
               ),
               column(10,
                      plotOutput("DEG_Violins_L2"),
                      dataTableOutput("DEG_L2_Table_Pair"))
             )),
    tabPanel("Cell-cell communication with CellChat",
             card(card_header("CellChat information flow plots"),
                  sidebarLayout(
                    sidebarPanel(
                      radioButtons("composite", h4("Select comparison to visualize for information flow plot"),
                                   choices = list("HC vs DCE", "HC vs DRE", "DCE vs DRE"),
                                   selected = "HC vs DCE"),
                      selectizeInput("pathways", h4("Select pathways to display"),
                                  choices = cellchatpathways, selected = "IL1",
                                  multiple = T),
                      actionButton("make_ranknet_plot", label = "Create plot")),
                    mainPanel(plotOutput("CellChatPlot"))
                  )),
             h4("CellChat circular plots for specific signaling pathways"),
             selectizeInput("pathway", h4("Select pathway to display"),
                            choices = cellchatpathways, selected = "IL1"),
             fluidRow(
               column(4, h4("HC"), hr(), plotOutput("CirclePlotHC")),
               column(4, h4("DCE"), hr(), plotOutput("CirclePlotDCE")),
               column(4, h4("DRE"), hr(), plotOutput("CirclePlotDRE"))),
             h4("CellChat ligand-receptor gene expression plot"),
             fluidRow(
               column(3,
                   radioButtons("comparison", h4("Select comparison to visualize"),
                                choices = list("HC vs DCE", "HC vs DRE", "DCE vs DRE"),
                                selected = "HC vs DCE"),
                   selectizeInput("lrpathway", h4("Select pathways to display"),
                               choices = cellchatpathways, selected = "IL1"),
                   actionButton("make_geneexp_plot", label = "Create plot")),
               column(9, plotOutput("cellchatgeneexpplot"))
               ))
  )
)
