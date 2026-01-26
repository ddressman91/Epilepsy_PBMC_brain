ui <- fluidPage(
  theme = shinythemes::shinytheme("united"),
  titlePanel("Brain sequencing data in epileptic patients and controls"),
  tabsetPanel(               
    tabPanel("Brain data overview",
             sidebarLayout(
               sidebarPanel(
                 radioButtons("BrainVar", h4("Select variable to visualize"),
                              choices = list("Cell Type","Region", "Diagnosis",
                                             "TCR clonotype group"),
                              selected = "Cell Type"),
                 downloadButton('downloadPlot_UMAP', 'Download Plot',size="s")
               ),
               mainPanel(plotOutput("UMAP")),
             )),
    tabPanel("Brain cell x gene",
             fluidRow(
               column(6,
                      radioButtons("CxG_Brain_UMAP",h3("Variable of interest"),
                                  choices= c("Cell Type","Region", "Diagnosis",
                                             "TCR clonotype group"),
                                  selected = "Cell Type"),
                      plotOutput("CxG_Brain_A")),
               column(6,
                      selectizeInput("CxG_marker_Brain",h3("Gene input"),
                                  choices = rownames(data_matrix),selected = "CD8A"),
                      plotOutput("CxG_FeaturePlot_Brain")))),
    tabPanel("Brain cluster marker genes",
             fluidRow(
               column(2,
                      selectizeInput("Vln_Brain_marker", h4("Cell type marker"),
                                  choices = rownames(data_matrix),selected = "CD8A"),
                      hr(),
                      selectInput("Vln_MG_celltype", h4("Microglial cluster for DGE table"),
                                     choices = paste0("MG",seq(1,9)),selected = "MG1"),
                      hr(),
                      downloadButton('downloadTable_VlnPlot_Brain', 'Download Table',size="s")),
               column(10,plotOutput("Vln_feature_Brain"))),
             hr(),
             dataTableOutput("DEG_MG_Table")),
    tabPanel("DRE vs control differential gene expression",
             h3("Select cell type"),
             fluidRow(
               column(2,
                      radioButtons("HC_TLE_Cluster",label="Cell type",
                                   choices=c("MG1", "MG2", "MG3", "MG4", "MG5", "MG6",
                                             "MG7","MG8","Monocytes","CD8T","NK")),
                      hr(),
                      actionButton("make_volcano_plot", label = "Create volcano plot"),
                      downloadButton('downloadPlots_HC_vs_TLE', 'Download violin plot',size="s"),
                      downloadButton('downloadTable_volc_results', 'Download table',size="s")
               ),
               column(10,girafeOutput("Volcano_HC_TLE")
               )),
             fluidRow(
               h4("Violin plots of top differentially expressed genes"),
               hr(),
               column(12,style = "height:450px",plotOutput("Vln_Plots_HC_TLE"))
             )),
    tabPanel("Cell-cell communication with CellChat",
             card(card_header("CellChat information flow plots"),
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("pathways", h4("Select pathways to display"),
                                     choices = brainpathways, selected = "IL1",
                                     multiple = T),
                      actionButton("make_ranknet_plot", label = "Create plot")),
                    mainPanel(plotOutput("BrainCellChatPlot"))
                  )),
             h4("CellChat circular plots for specific signaling pathways"),
             selectizeInput("pathway", h4("Select pathway to display"),
                            choices = brainpathways, selected = "IL1"),
             fluidRow(
               column(6, h4("Control"), hr(), plotOutput("BrainCirclePlotCtrl")),
               column(6, h4("DRE"), hr(), plotOutput("BrainCirclePlotTLE"))),
             h4("CellChat ligand-receptor gene expression plot"),
             fluidRow(
               column(3,
                      selectizeInput("brainlrpathway", h4("Select pathways to display"),
                                     choices = brainpathways, selected = "IL1"),
                      actionButton("make_geneexp_plot", label = "Create plot")),
               column(9, plotOutput("braincellchatgeneexpplot"))
             ))
    ))
