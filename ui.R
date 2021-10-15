library(shiny)
library(plotly)
library(ggplot2)
library(shinyjs)
library(DT)
library(Seurat)
library(SeuratData)
library(cowplot)
library(sctransform)
#library(SeuratWrappers)
library(dplyr)
library(pbmcapply)
#library(harmony)
#library(rliger)
library(reshape2)
library(shinydashboard)
library(shinyalert)
library(shinyFiles)
library(shinyWidgets)

shiny_one_panel = fluidPage(
    titlePanel("Seurat analysis of scRNAseq data"),
    hr(),

    fluidRow(
        ##------Sidebar---------
        column(3,
               h4('Load Data:'),
               wellPanel(
                   fileInput(inputId = 'tpmFiles',
                             label = "Gene expression file",
                             multiple = FALSE,
                             accept = c("text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".Robj")),
                   fileInput(inputId = 'cellAnnoFiles',
                             label = "Metadata",
                             multiple = FALSE,
                             accept = c("text/csv",
                                        "text/comma-separated-values,text/plain")),
                   checkboxInput(inputId = 'norm',
                                 label = "Normalise?",
                                 value = TRUE),
                   fluidRow(
                       column(6,
                              numericInput(inputId = "min.genes",
                                           label = "Min. genes",
                                           value = 200,
                                           min = 1)
                       ),
                       column(6,
                              numericInput(inputId = "min.cells",
                                           label = "Min. cells",
                                           value = 3,
                                           min = 1)
                       )
                   ),
                   textInput(inputId = "projName",
                             label = "Project Name",
                             value = "Seurat_analysis"),
                   fluidRow(
                       column(6,
                              actionButton("loadButton", "Create  Seurat Object", icon = icon("hand-o-right"))
                       ),
                       column(6,
                              actionButton("reset", "Reset Data", icon = icon("repeat"))
                       )
                   )
               ),

               ##------Plot download---------
               h4("Export to PDF:"),
               wellPanel(
                   ## Conditional panel for different plots
                   conditionalPanel(" input.QC == 'QC_panel1' && input.tabs == 'QC plots' ",
                                    actionButton("PDFa", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.QC == 'QC_panel2' && input.tabs == 'QC plots' ",
                                    actionButton("PDFb", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.tabs == 'Normalization and Variable Gene Plot' ",
                                    actionButton("PDFc", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.Pca == 'P_panel1' && input.tabs == 'PCA' ",
                                    actionButton("PDFd", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.Pca == 'P_panel2' && input.tabs == 'PCA' ",
                                    actionButton("PDFe", "Download", icon = icon("download"))
                   ),
                   #conditionalPanel(" input.Sct == 'SCT_panel' && input.tabs == 'SCT' ",
                   #                    actionButton("PDFs", "Download", icon = icon("download"))
                   #),
                   conditionalPanel(" input.Pca == 'P_panel3' && input.tabs == 'PCA' ",
                                    actionButton("PDFg", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.Pca == 'P_panel4' && input.tabs == 'PCA' ",
                                    actionButton("PDFh", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.tabs == 'Clustering' ",
                                    actionButton("PDFf", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.tabs == 'UMAP' ",
                                    actionButton("PDFi", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.tabs == 'TSNE' ",
                                    actionButton("PDFj", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.tabs == 'DEGs' ",
                                    actionButton("PDFk", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.tabs == 'Data Integration using Seurat' ",
                                    actionButton("PDFl", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.tabs == 'Data Integration using Harmony' ",
                                    actionButton("PDFm", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.tabs == 'Data Integration using RLiger' ",
                                    actionButton("PDFn", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.tabs == 'scMultiomics using Seurat' ",
                                    actionButton("PDFo", "Download", icon = icon("download"))
                   ),
                   ## ensure no spill over in button text
                   tags$head(
                       tags$style(HTML('
                                   .btn {
                                   white-space: normal;
                                   }'
                       )
                       )
                   ),
                   ## Conditional is separate from pdf options
                   hr(),
                   fluidRow(
                       column(6,
                              sliderInput(inputId="pdf_w", label = "PDF width(in):",
                                          min=3, max=20, value=8, width=100, ticks=F)
                       ),
                       column(6,
                              sliderInput(inputId="pdf_h", label = "PDF height(in):",
                                          min=3, max=20, value=8, width=100, ticks=F)
                       )),

                   #actionButton("OpenDir", "Open download folder", icon = icon("folder"))
               ),

               ##------Save Data---------
               hr(),
               actionButton("saveButton", "Save Data", icon = icon("hand-o-right")),

               hr(),
               h4(tags$a(href="mailto:Chen_Jinmiao@immunol.a-star.edu.sg?subject=[cytof-question]",
                         "Contact Us")),
               imageOutput("logo", height = "60px")
        ),
        ##------Main area---------
        column(9,
               tabsetPanel(type = "pills", id = "tabs",
                           ## add file preview tab
                           ##------QC plots---------
                           tabPanel("QC plots", fluidPage(
                               hr(),
                               tabsetPanel(id="QC",
                                           tabPanel(title="Violin Plots", value = "QC_panel1",
                                                    br(),
                                                    fluidRow(
                                                        column(6,
                                                               plotlyOutput("nFeature_RNAPlot", width = "100%"),
                                                               br(),
                                                               plotlyOutput("mitoPlot", width = "100%"),
                                                               br(),
                                                               plotlyOutput("nCount_RNAPlot", width = "100%")
                                                        ),
                                                        #column(6,
                                                        #       verbatimTextOutput("name")
                                                        #)
                                                    )

                                           ),
                                           tabPanel(title="Feature Scatter Plots", value="QC_panel2",
                                                    br(),
                                                    fluidRow(
                                                        column(6,
                                                               plotlyOutput("FeatureScatterPlot1", width = "100%"),
                                                               br(),
                                                               plotlyOutput("FeatureScatterPlot2", width = "100%")
                                                        )
                                                    )
                                           )
                               )
                           )),
                           ##------Normalization and Variable Genes---------
                           tabPanel("Normalization and Variable Gene Plot", fluidPage(
                               hr(),

                               selectInput("norm1",
                                           label = "Normalization method",
                                           choices = c("LogNormalize", "CLR", "RC")
                               ),

                               textOutput("nVarGenes"),
                               fluidRow(
                                   column(4,
                                          numericInput("y.cutoff",
                                                       label = "Standardized Variance",
                                                       value = 0.5)
                                   ),
                                   column(4,
                                          numericInput("x.cutoff",
                                                       label = "Average Expression",
                                                       value = 0.1,
                                                       min = 0)
                                   ),
                                   column(4,
                                          actionButton("findVarGenes", "Analyse variable genes", icon = icon("hand-pointer-o")),
                                          actionButton("doSCTransform", "Run SCTransform", icon = icon("hand-pointer-o"))
                                          # actionButton("doVarplot", "Plot variable genes", icon = icon("hand-pointer-o"))
                                   )),
                               plotOutput("VarGenes", width = "100%")
                           )),

                           ##------PCA---------
                           tabPanel("PCA", fluidPage(
                               hr(),
                               tabsetPanel(id="Pca",
                                           tabPanel(title="PCA Plot", value="P_panel1",
                                                    br(),
                                                    fluidRow(
                                                        column(3,
                                                               actionButton("doPCA", "Run PCA", icon = icon("hand-pointer-o"))
                                                        ),
                                                        #   column(3,
                                                        #          selectInput("x.pc",
                                                        #                      label = "X-axis PC to use",
                                                        #                      choices = 1:20,
                                                        #                      selected = 1)
                                                        #placeholder for z-axis input
                                                        #   ),
                                                        #   column(3,
                                                        #          selectInput("y.pc",
                                                        #                      label = "Y-axis PC to use",
                                                        #                      choices = 1:20,
                                                        #                      selected = 2)
                                                        #   ),
                                                        #   column(3,
                                                        #          selectInput("z.pc",
                                                        #                      label = "Z-axis PC to use",
                                                        #                      choices = 1:20,
                                                        #                      selected = 3)
                                                        #   )
                                                        # ),
                                                        # fluidRow(
                                                        #     column(3 ,
                                                        #            numericInput("pc.plot.size",
                                                        #                         label = "Point Size:",
                                                        #                         value = 1,
                                                        #                         min = 0.1,
                                                        #                         step = 0.5)
                                                        #     ),
                                                        #     column(3,
                                                        #            sliderInput("pca.plot.alpha",
                                                        #                        label = "Point Transparency",
                                                        #                        min = 0,
                                                        #                        max = 1,
                                                        #                        step = 0.1,
                                                        #                        value = 0.8)
                                                        #     ),
                                                        # column(6,
                                                        #       uiOutput("clustUI")
                                                        #)
                                                    ),
                                                    selectInput("assays1",
                                                                label = "Normalize by:",
                                                                choices = c("RNA", "SCT")
                                                    ),
                                                    plotlyOutput("PCA2DPlot", width = "100%"),
                                                    #plotlyOutput("PCA3DPlot", width = "100%")
                                           ),
                                           tabPanel(title="PC Gene Visualisation", value="P_panel2",
                                                    br(),
                                                    selectInput("select.pc",
                                                                label = "PC to plot",
                                                                choices = c(1:20)
                                                    ),
                                                    fluidRow(
                                                        column(4,
                                                               plotOutput("vizPlot", width = "100%", height = "600px")
                                                        ),
                                                        column(8,
                                                               plotOutput("PCHeatmap", width = "100%", height = "600px")
                                                        )
                                                    ),
                                                    DT::dataTableOutput("PCtable")
                                           ),
                                           tabPanel(title="JackStraw", value="P_panel3",
                                                    br(),
                                                    actionButton("doJack", label = "Run Jackstraw"),
                                                    br(),
                                                    plotOutput("Jackstraw", width = "100%")
                                           ),
                                           tabPanel(title="Elbow", value="P_panel4",
                                                    br(),
                                                    actionButton("doElbow", label = "Get Elbow Plot"),
                                                    br(),
                                                    br(),
                                                    plotOutput("Elbow", width = "100%")

                                           )
                               )
                           )),
                           tabPanel("Clustering", fluidPage(
                               hr(),
                               fluidRow(
                                   column(3,
                                          numericInput("clus.res",
                                                       label = "Resolution used",
                                                       value = 0.6,
                                                       min = 0.1,
                                                       step = 0.1)
                                   ),
                                   br(),
                                   selectInput("dim.used",
                                               label = "PC to plot",
                                               choices = c(2:50)
                                   ),
                                   column(3,
                                          br(),
                                          actionButton("findCluster", "Find Clusters", icon = icon("hand-pointer-o")),
                                          textOutput("cluster.done"),
                                          br()
                                   )),
                               br(),
                               plotlyOutput("Cluster2DPlot_1", width = "100%"),
                               br(),
                               plotlyOutput("Cluster2DPlot_2", width = "100%"),
                               br(),
                               plotlyOutput("Cluster2DPlot_3", width = "100%"),
                               #plotlyOutput("Umap_3d_plot", width = "100%"),
                               fluidRow(
                               )
                           )),



                           ##------UMAP---------
                           tabPanel("UMAP", fluidPage(
                               hr(),
                               fluidRow(
                                   column(3,
                                          numericInput("dim.used",
                                                       label = "Dimensions used",
                                                       value = 10)
                                   ),
                                   br(),
                                   #column(3,
                                   #      numericInput("max.iter1",
                                   #                   label = "Max Iterations",
                                   #                   value = 2000,
                                   #                   min = 100)
                                   #),
                                   column(3,
                                          br(),
                                          actionButton("doUmap", "Run UMAP", icon = icon("hand-pointer-o")),
                                          textOutput("Umap.done"),
                                          br()
                                   )),
                               #  fluidRow(
                               #      column(3,
                               #             numericInput("umap.plot.size",
                               #                          label = "Point Size:",
                               #                          value = 1,
                               #                          min = 0.1,
                               #                          step = 0.1)
                               #      ),
                               #      column(3,
                               #             sliderInput("umap.plot.alpha",
                               #                         label = "Point Transparency",
                               #                         min = 0,
                               #                         max = 1,
                               #                         step = 0.1,
                               #                         value = 0.8)
                               #      ),
                               #      column(2,
                               #             textOutput("selection.summary1"),
                               #             textInput("selection.name", label = "New cluster name", value = "custom")
                               #      ),
                               #column(4,
                               #      actionButton("create.selection1", label = "Create cluster from selection"),
                               #      actionButton("reset.selection1", label = "Reset identities")
                               #   )
                               #),
                               br(),
                               plotlyOutput("Umap_2d_plot_1", width = "100%"),
                               br(),
                               plotlyOutput("Umap_2d_plot_2", width = "100%"),
                               br(),
                               plotlyOutput("Umap_2d_plot_3", width = "100%"),
                               #plotlyOutput("Umap_3d_plot", width = "100%"),
                               fluidRow(
                               )
                           )),
                           ##------TSNE---------
                           tabPanel("TSNE", fluidPage(
                               hr(),
                               fluidRow(
                                   column(3,
                                          numericInput("dim.used",
                                                       label = "Dimensions used",
                                                       value = 10)
                                   ),
                                   br(),
                                   #column(3,
                                   #          numericInput("max.iter",
                                   #                       label = "Max Iterations",
                                   #                       value = 2000,
                                   #                       min = 100)
                                   #   ),
                                   column(3,
                                          uiOutput("perplex.option")
                                   ),
                                   column(3,
                                          br(),
                                          actionButton("doTsne", "Run TSNE", icon = icon("hand-pointer-o")),
                                          textOutput("Tsne.done"),
                                          br()
                                   )),
                               #fluidRow(
                               #   column(3,
                               #          numericInput("tsne.plot.size",
                               #                       label = "Point Size:",
                               #                       value = 1,
                               #                       min = 0.1,
                               #                       step = 0.1)
                               #   ),
                               #   column(3,
                               #          sliderInput("tsne.plot.alpha",
                               #                      label = "Point Transparency",
                               #                      min = 0,
                               #                      max = 1,
                               #                      step = 0.1,
                               #                      value = 0.8)
                               #   ),
                               #   column(2,
                               #          textOutput("selection.summary"),
                               #          textInput("selection.name", label = "New cluster name", value = "custom")
                               #   ),
                               #   column(4,
                               #          actionButton("create.selection", label = "Create cluster from selection"),
                               #          actionButton("reset.selection", label = "Reset identities")
                               #   )
                               # ),
                               br(),
                               plotlyOutput("Tsne_2d_plot_1", width = "100%"),
                               br(),
                               plotlyOutput("Tsne_2d_plot_2", width = "100%"),
                               br(),
                               plotlyOutput("Tsne_2d_plot_3", width = "100%"),
                               #plotlyOutput("Tsne_3d_plot", width = "100%"),
                               fluidRow(
                               )
                           )),
                           ##------DEGs---------
                           tabPanel("DEGs", fluidPage(
                               hr(),
                               fluidRow(
                                   #   column(4,
                                   #          uiOutput("clust1")
                                   #   ),
                                   #   column(4,
                                   #          uiOutput("clust2")
                                   #   ),

                                   column(3, selectInput("min_pct",
                                                         label = "min.pct",
                                                         choices = c("0.1", "0.25"))
                                   ),

                                   column(3, selectInput("logfc",
                                                         label = "logfc.threshold",
                                                         choices = c("0.1", "0.25"))
                                   ),

                                   column(4,
                                          actionButton("doDeg", "Run DEGs", icon = icon("hand-pointer-o"))
                                   )),
                               fluidRow(
                                   column(6,
                                          uiOutput("deg.gene.select"),
                                          plotlyOutput("Deg.plot", width = "100%"),
                                          br(),
                                          plotlyOutput("Deg1.plot", width = "100%"),
                                          br(),
                                          plotlyOutput("Deg2.plot", width = "100%")

                                   ),
                                   column(6,
                                          DT::dataTableOutput("Deg.table"),
                                          br(),
                                          plotlyOutput("Deg3.plot", width = "100%")
                                   ))
                           ))




                           ##------END---------
               )
        )
    )
)

dbHeader <- dashboardHeader(title = "cytofkit2")
dbHeader$children[[2]]$children <-  tags$a(href='https://github.com/JinmiaoChenLab',
                                           tags$img(src='https://avatars1.githubusercontent.com/u/8896007?s=400&u=b0029c2e64f405ea0a46d311239b674a430ec77c&v=4'
                                                    ,height='60',width='60', align='left')
                                           , tags$div('cytofkit2', style='color:white;font-family:arial rounded MT bold'))

dashboardPage(skin = "yellow",
              dbHeader,
              dashboardSidebar(
                  sidebarMenu(id = "sbm",
                              menuItem(tags$p(style = "display:inline;font-size: 20px;", "Seurat"), tabName = "seurat", icon = icon('cog'))


                  )# end of sidebarMenu
              ),#end of dashboardSidebar
              dashboardBody(
                  #includeCSS("www/custom.css")
                  useShinyalert()
                  , shinyjs::useShinyjs()
                  , tabItem(
                      tabName = "seurat"
                      , shiny_one_panel
                  ) # End of tabItem

              )# end of dashboard body
)# end of dashboard page
