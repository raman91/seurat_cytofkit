## max data size
options(shiny.maxRequestSize = 1024^10)
options(shiny.launch.browser = T)

shinyServer(function(input, output, session) {
  v <- reactiveValues(scData = NULL,
                      scData1 = NULL,
                      scData2 = NULL,
                      scData3 = NULL,
                      scData4 = NULL,
                      scData5 = NULL,
                      scData6 = NULL,
                      scData7 = NULL,
                      scData8 = NULL,
                      scData9 = NULL,
                      scData10 = NULL,
                      scData11 = NULL,
                      scData12 = NULL,
                      scDatat = NULL,
                      idents = NULL,
                      isPCAdone = NULL,
                      isUMAPdone = NULL,
                      isTSNEdone = NULL,
                      isPCAdone1 = NULL,
                      isUMAPdone1 = NULL,
                      isTSNEdone1 = NULL,
                      isPCAdone2 = NULL,
                      isUMAPdone2 = NULL,
                      isTSNEdone2 = NULL,
                      isPCAdone3 = NULL,
                      isUMAPdone3 = NULL,
                      isTSNEdone3 = NULL,
                      isUMAPdone4 = NULL,
                      isTSNEdone4 = NULL,
                      isVisdone = NULL,
                      isClusterdone = NULL,
                      isDataIntegration = NULL,
                      pcGenes = NULL,
                      plotlySelection = NULL,
                      ips.markers = NULL)
  #celltypes <- NULL
  prePlot <- function(){
    while(names(dev.cur()) != "null device"){
      dev.off()
    }
  }
  observe({
    #s <- event_data("plotly_selected")
    #cells <- s[["key"]]
    v$plotlySelection <- event_data("plotly_selected")[["key"]]
  })
  ##-------------------Side Panel-------------------

  normMethod <- NULL

  output$name.field <- renderUI({
    if(is.null(input$cellAnnoFiles)){
      numericInput(inputId = "field",
                   label = "Field",
                   value = 1,
                   min = 1)
    }else{
      annoFile <- input$cellAnnoFiles
      anno.data <- read.table(annoFile$datapath[1], header = T,
                              sep = "\t", stringsAsFactors = FALSE)
      groupings <- colnames(anno.data)
      selectInput("groupby",
                  label = "Group by:",
                  choices = groupings)
    }
  })

  observeEvent(input$loadButton, {
    tpmFiles <- input$tpmFiles
    annoFile <- input$cellAnnoFiles
    names.field <- input$field
    if (is.null(tpmFiles)){
      v$scData <- NULL
    }else{
      withProgress(message="Loading and Processing Data...", value=0, {
        print(tpmFiles$datapath)
        print(tpmFiles$name)
        print(file.exists(paste(tpmFiles$datapath[1], "/", tpmFiles$name[1], sep="")))
        exp.data <- read.table(tpmFiles$datapath,
                               sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)
        #additional.ident1 <- NULL
        if(!is.null(annoFile)){
          anno.data <- read.table(annoFile$datapath[1], header = T,
                                  sep = "\t", stringsAsFactors = FALSE, row.names=1)
          #to.append <- apply(anno.data1, 1, paste, collapse = "_")
          #colnames(exp.data1) <- to.append
          #names.field1 <- match(input$groupby, colnames(anno.data1))
          #additional.ident1 <- data.frame(data.frame(anno.data1[,-1], row.names = to.append))
          #additional.ident1[] <- lapply(additional.ident1, factor)
        }
        incProgress(0.5, "Creating Seurat Object")
        sObj <- CreateSeuratObject(exp.data,
                                   meta.data = anno.data,
                                   project = input$projName,
                                   names.field = names.field,
                                   names.delim = input$delim,
                                   is.expr = input$expThres,
                                   normalization.method = normMethod,
                                   min.genes = input$min.genes,
                                   min.cells = input$min.cells)
        #mito.genes <- grep("^MT-", rownames(sObj@assays$RNA@data), ignore.case = TRUE, value = TRUE)
        sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
        #incProgress(0.5, "Adding metadata")
        #sObj <- AddMetaData(sObj, percent.mt, "percent.mt")
        #if(!is.null(additional.ident1)){
        #  sObj1 <- AddMetaData(sObj1, additional.ident1)
        #}
        v$scData <- sObj
      })
    }
  })
  dir.create("Seurat_results")
  #})

  observeEvent(input$reset, {
    session$reload()
    print("Reset done")
  })

  observeEvent(input$saveButton, {
    if(!is.null(input$tpmFiles)){
      withProgress(message="Saving Results...", value=0, {
        print(getwd())
        dir.create("Seurat_results")
        resultDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results")
        filename <- paste0(resultDir, .Platform$file.sep, v$scData@project.name, "_", Sys.Date())
        sObj <- v$scData
        save(sObj, file= paste0(resultDir, .Platform$file.sep, sObj@project.name, "_", Sys.Date(), ".Robj"))
      })
      ## open the results directory
      opendir(resultDir)
    }
  })

  output$ident.swap <- renderUI({
    if(is.null(v$scData)){
      return(NULL)
    }else{
      groupings1 <- names(v$scData@meta.data[,!names(v$scData@meta.data) %in% c("nFeature_RNA", "nCount_RNA", "percent.mt")])
      tagList(
        h4("Set current identity:"),
        fluidRow(
          column(6,
                 selectInput("active.ident", label = NULL,
                             choices = groupings1)
          ),
          column(6,
                 actionButton("swap.ident",label = NULL, icon = icon("arrow-right"))
          )
        )

      )
    }
  })

  observeEvent(input$swap.ident, {
    v$scData <- SetIdent(v$scData, value = as.character(v$scData@meta.data[,input$active.ident]))
  })

  output$logo <- renderImage({
    return(list(
      src = "inst/extdata/logo.png",
      contentType = "image/png",
      alt = "Singapore Immunology Network"
    ))
  }, deleteFile = FALSE)

  opendir <- function(dir = getwd()){
    if (.Platform['OS.type'] == "windows"){
      shell.exec(dir)
    } else {
      system(paste(Sys.getenv("R_BROWSER"), dir))
    }
  }

  observeEvent(input$OpenDir, {
    resultDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results")
    if(!dir.exists(resultDir)){
      dir.create("Seurat_results")
    }
    pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
    if(dir.exists(pdfDir)){
      opendir(pdfDir)
    }else{
      warning("No reports created yet!")
      dir.create(pdfDir)
    }
  })

  ##---------------QC tabset-------------------

  output$nFeature_RNAPlot <- renderPlotly({
    if(is.null(v$scData)){
      plotly_empty()
    }else{
      VlnPlot(v$scData, "nFeature_RNA")
    }
  })

  output$mitoPlot <- renderPlotly({
    if(is.null(v$scData)){
      plotly_empty()
    }else{
      VlnPlot(v$scData, "percent.mt")
    }
  })

  output$nCount_RNAPlot <- renderPlotly({
    if(is.null(v$scData)){
      plotly_empty()
    }else{
      VlnPlot(v$scData, "nCount_RNA")
    }
  })

  output$name <- renderPrint({
    s <- event_data("plotly_selected")
    c(s[["key"]], class(s[["key"]]))
  })

  observeEvent(input$PDFa, {
    if(!is.null(v$scData)){
      withProgress(message="Downloading plot PDF files...", value=0, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"QC_violin_plot_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,
                              "QC_violin_plot_",
                              Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        prePlot()
        nG <- VlnPlot(v$scData, "nFeature_RNA")
        pM <- VlnPlot(v$scData, "percent.mt")
        nU <- VlnPlot(v$scData, "nCount_RNA")
        pdf(filename2,
            width=as.numeric(input$pdf_w),
            height=as.numeric(input$pdf_h))
        print(nG)
        print(pM)
        print(nU)
        dev.off()
      })
    }
  })

  ## FeatureScatter plot

  output$FeatureScatterPlot1 <- renderPlotly({
    if(is.null(v$scData)){
      plotly_empty()
    }else{
      print(FeatureScatter(v$scData, "nCount_RNA", "nFeature_RNA"))
    }
  })

  output$FeatureScatterPlot2 <- renderPlotly({
    if(is.null(v$scData)){
      plotly_empty()
    }else{
      print(FeatureScatter(v$scData, "nCount_RNA", "percent.mt"))
    }
  })

  observeEvent(input$PDFb, {
    if(!is.null(v$scData)){
      withProgress(message="Downloading plot PDF files...", value=0, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"QC_featurescatter_plot_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,
                              "QC_featurescatter_plot_",
                              Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        prePlot()
        pdf(filename2,
            width=as.numeric(input$pdf_w),
            height=as.numeric(input$pdf_h))
        feature1 <- FeatureScatter(v$scData, "nCount_RNA", "nFeature_RNA")
        feature2 <- FeatureScatter(v$scData, "nCount_RNA", "percent.mt")
        print(feature1)
        print(feature2)
        dev.off()
      })
    }
  })

  ##---------------Normalization and Variable Genes tabset-------------------

  observeEvent(input$findVarGenes, {
    withProgress(message = "Finding variable genes...", value = 0, {
      v$scData <- NormalizeData(v$scData, normalization.method = input$norm1)
      v$scData <- FindVariableFeatures(v$scData,
                                       mean.function = ExpMean,
                                       dispersion.function = LogVMR)
      #all.genes <- rownames(v$scData)
      v$scData <- ScaleData(v$scData)
      incProgress(0.5)
      VarGeneText <- paste0("Number of variable genes: ", length(v$scData@assays$RNA@var.features))
      output$nVarGenes <- renderText(VarGeneText)
      varGenePlotInput <- function(){
        if(is.null(v$scData)){
          return(NULL)
        }else{
          withProgress(message="Plotting variable genes...", value=0, {
            top10 <- head(VariableFeatures(v$scData), 10)
            variable_feature1 <- VariableFeaturePlot(v$scData)
            variable_feature2 <- LabelPoints(plot = variable_feature1, points = top10, repel = TRUE)
            print (variable_feature1)
            print (variable_feature2)
            dev.off()
          })
          #  observeEvent(input$doVarplot,{

          #  })
        }
      }

      observeEvent(input$doSCTransform, {
        withProgress(message = "Performing scTransform...", value = 0,{
          incProgress(0.5, message = "Running scTransform...")
          v$scData <- SCTransform(v$scData, vars.to.regress = "percent.mt", verbose = FALSE)
        })
      })

      output$VarGenes <- renderPlot({
        varGenePlotInput()
      }, height = 800, width = 850)
      observeEvent(input$PDFc, {
        if(!is.null(v$scData)){
          withProgress(message="Downloading plot PDF files...", value=0, {
            print(getwd())
            pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
            if(!dir.exists(pdfDir)){
              dir.create(pdfDir)
            }
            filename2 <- paste0(pdfDir, .Platform$file.sep,"Var_genes_plot_", Sys.Date(), ".pdf")
            i = 0
            while(file.exists(filename2)){
              filename2 <- paste0(pdfDir, .Platform$file.sep,
                                  "Var_genes_plot_",
                                  Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
              i = i + 1;
            }
            prePlot()
            pdf(filename2,
                width=as.numeric(input$pdf_w),
                height=as.numeric(input$pdf_h))
            plot1 <- VariableFeaturePlot(v$scData)
            print(plot1)
            dev.off()
            txtfile <- sub("Var_genes_plot_", "Var_gene_list_", filename2)
            txtfile <- sub(".pdf", ".txt", txtfile)
            write(v$scData@assays$RNA@var.features, file = txtfile)
          })
        }
      })
    })
  })

  ##---------------PCA tabset-------------------
  # PCA plot
  observeEvent(input$doPCA, {
    withProgress(message = "Scaling Data...", value = 0,{
      incProgress(0.5, message = "Running PCA...")
      v$scData <- RunPCA(v$scData, features = VariableFeatures(object = v$scData), assay = input$assays1)
      print(v$scData[["pca"]], dims = 1:5, nfeatures = 5)
      v$isPCAdone <- TRUE
      PCA_plot <- DimPlot(v$scData, reduction = "pca", label = T)
      print(PCA_plot)
      incProgress(0.4, message = "Getting list of PC genes...")
      pc.table <- list()
      for(i in 1:20){
        pcg <- TopCells(v$scData)
        pc.table[[i]] <- pcg
      }
      pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
      v$pcGenes <- pc.table
    })
  })

  output$PCA2DPlot <- renderPlotly({
    if(is.null(v$isPCAdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating 2D PCA Plot...", value=0, {
        DimPlot(v$scData, reduction = "pca", label = T)
      })
    }
  })

  #output$PCA3DPlot <- renderPlotly({
  #    if(is.null(v$isPCAdone)){
  #        plotly_empty()
  #    }else{
  #        withProgress(message="Generating 3D PCA Plot...", value=0, {
  #            DimPlot(v$scData, reduction = "pca", label = T)
  #        })
  #    }
  #})

  observeEvent(input$PDFd, {
    if(!is.null(v$scData)){
      withProgress(message="Downloading plot PDF files...", value=0, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"PCA_plot_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,
                              "PCA_plot_",
                              Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        pcaplot <- DimPlot(v$scData, reduction = "pca", label = T)
        prePlot()
        pdf(filename2,
            width=as.numeric(input$pdf_w),
            height=as.numeric(input$pdf_h))
        print(pcaplot)
        dev.off()
      })
      withProgress(message="Downloading PCA coordinates...", value=0.5, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"pca_", Sys.Date(), ".txt")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,
                              "pca_",
                              Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        write.csv(v$scData@reductions$pca@cell.embeddings, file = filename2)
      })
      withProgress(message="Downloading cluster IDs...", value=0.9, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"cluster_", Sys.Date(), ".txt")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,"cluster_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        write.csv(v$scData@active.ident, file = filename2)
      })
    }
  })

  # Viz plot

  output$vizPlot <- renderPlot({
    if(is.null(v$scData)){
      return(NULL)
    }else{
      VizDimLoadings(v$scData, dims = as.numeric(input$select.pc))
    }
  })

  output$PCHeatmap <- renderPlot({
    if(is.null(v$scData)){
      return(NULL)
    }else{
      DimHeatmap(v$scData, dims = as.numeric(input$select.pc))
    }
  })

  output$PCtable <- DT::renderDataTable({
    if(is.null(v$scData) ){
      return(NULL)
    }else{
      v$pcGenes
    }
  }, options = list(scrollX = TRUE))

  observeEvent(input$PDFe, {
    if(!is.null(v$scData)){
      withProgress(message="Downloading plot PDF files...", value=0, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"Viz_Heatmap_plots_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,
                              "Viz_Heatmap_plots_",
                              Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        prePlot()
        pdf(filename2,
            width=as.numeric(input$pdf_w),
            height=as.numeric(input$pdf_h))
        #isolate({
        vizdim <- VizDimLoadings(v$scData, pcs.use = as.numeric(input$select.pc))
        dimheat <- DimHeatmap(v$scData)
        print (vizdim)
        print (dimheat)
        #})
        dev.off()
        pcGenes <- v$pcGenes
        write.csv(v$pcGenes, file = paste0(pdfDir, .Platform$file.sep,"PC_genes_", Sys.Date(), ".csv"))
      })
    }
  })

  # Jackstraw
  observeEvent(input$doJack, {
    if(is.null(v$scData)){
      return(NULL)
    }else{
      withProgress(message="Running Jackstraw...", value=0.5, {
        v$scData <- JackStraw(v$scData, num.replicate = 100)
        v$scData <- ScoreJackStraw(v$scData, dims = 1:20)
        #selected.pc <- which(v$scData@reductions$pca@jackstraw@overall.p.values < 0.01)
        #v$signPC <- selected.pc
      })
    }
  })

  output$Jackstraw <- renderPlot({
    if(is.null(v$scData)){
      return(NULL)
    }else{
      JackStrawPlot(v$scData, dims = 1:20)
    }
  }, height = 800, width = 850)

  observeEvent(input$PDFg, {
    if(!is.null(v$scData)){
      withProgress(message="Downloading plot PDF files...", value=0, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"Jackstraw_plot_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,"Jackstraw_plot_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        prePlot()
        pdf(filename2,
            width=as.numeric(input$pdf_w),
            height=as.numeric(input$pdf_h))
        jack <- JackStrawPlot(v$scData, dims = 1:20)
        print (jack)
        dev.off()
      })
    }
  })

  # Elbow
  output$Elbow <- renderPlot({
    if(is.null(v$scData@reductions$pca@jackstraw)){
      return(NULL)
    }else{
      withProgress(message="Generating Elbow Plot...", value=0.5, {
        ElbowPlot(v$scData)
      })
    }
  }, height = 800, width = 850)

  observeEvent(input$PDFh, {
    if(!is.null(v$scData)){
      withProgress(message="Downloading plot PDF files...", value=0, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"Elbow_plot_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,"Elbow_plot_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        prePlot()
        pdf(filename2,
            width=as.numeric(input$pdf_w),
            height=as.numeric(input$pdf_h))
        print(ElbowPlot(v$scData))
        dev.off()
      })
    }
  })

  ##---------------Clustering-------------------

  output$clustUI <- renderUI({
    if(is.null(v$isPCAdone)){
      return(NULL)
    }else{
      tagList(
        fluidRow(
          column(6,
                 numericInput("clus.res",
                              label = "Cluster Resolution",
                              value = 0.6,
                              min = 0.1,
                              step = 0.1)
          ),
          br(),
          selectInput("dim.used",
                      label = "First n PCs",
                      choices = c(2:50)
          ),
          column(6,
                 actionButton("findCluster", "Find Clusters", icon = icon("hand-pointer-o")),
                 textOutput("cluster.done")
          )
        )
      )
    }
  })

  observeEvent(input$findCluster, {
    withProgress(message = "Finding clusters...", value = 0.3, {
      v$scData <- FindNeighbors(v$scData, dims = 1:input$dim.used, assay = input$assays1)
      v$scData <- FindClusters(v$scData, resolution = input$clus.res)
      output$cluster.done <- renderText(paste0("Clustering done!"))
      v$isClusterdone <- TRUE
    })
  })

  output$Cluster2DPlot_1 <- renderPlotly({
    if(is.null(v$isClusterdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating 2D Cluster Plot...", value=0, {
        DimPlot(v$scData, reduction = "pca", label = T)
      })
    }
  })

  output$Cluster2DPlot_2 <- renderPlotly({
    if(is.null(v$isClusterdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating 2D Cluster Plot...", value=0, {
        DimPlot(v$scData, reduction = "pca", label = T, group.by = 'batch')
      })
    }
  })

  output$Cluster2DPlot_3 <- renderPlotly({
    if(is.null(v$isClusterdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating 2D Cluster Plot...", value=0, {
        DimPlot(v$scData, reduction = "pca", label = T, group.by = 'celltype')
      })
    }
  })

  ##---------------UMAP tabset-------------------

  observeEvent(input$doUmap, {
    withProgress(message = "Running UMAP...", value = 0.3, {
      #dims.use <- NULL
      #if(!is.null(v$scData@reductions$pca@jackstraw)){
      #score.df <- JackStraw_pval(v$scData)
      #dims.use <- signif_PCs(score.df)
      #}else{
      #dims.use <- 1:input$dim.used
      #}
      v$scData <- RunUMAP(v$scData, dims = 1:input$dim.used, assay = input$assays1)
      output$Umap.done <- renderText(paste0("UMAP done!"))
      v$isUMAPdone <- TRUE
    })
  })

  output$Umap_2d_plot_1 <- renderPlotly({
    if(is.null(v$scData) || is.null(v$isUMAPdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
        DimPlot(v$scData, reduction = "umap", label = T)
      })
    }
  })

  output$Umap_2d_plot_2 <- renderPlotly({
    if(is.null(v$scData) || is.null(v$isUMAPdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
        DimPlot(v$scData, reduction = "umap", label = T, group.by = 'batch')
      })
    }
  })

  output$Umap_2d_plot_3 <- renderPlotly({
    if(is.null(v$scData) || is.null(v$isUMAPdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
        DimPlot(v$scData, reduction = "umap", label = T, group.by = 'celltype')
      })
    }
  })

  # output$Umap_3d_plot <- renderPlotly({
  #   if(is.null(v$scData) || is.null(v$isUMAPdone)){
  #     plotly_empty()
  #  }else{
  #     withProgress(message="Generating UMAP 3D Plot...", value=0, {
  #     DimPlot(v$scData, reduction = "umap", label = T)
  #      })
  #    }
  # })

  output$selection.summary1 <- renderText({
    if(is.null(v$plotlySelection)){
      return(NULL)
    }else{
      t1 <- paste0(length(v$plotlySelection), " cells selected")
      t1
    }
  })

  observeEvent(input$create.selection1, {
    ## stash old identity
    if(is.null(v$scData@meta.data$cluster.ident)){
      v$scData <- StashIdent(object = v$scData, save.name = 'cluster.ident')
    }
    v$scData <- Idents (object = v$scData,
                        cells = v$plotlySelection,
                        value = as.character(input$selection.name)
    )
    updateTabsetPanel(session, "tabs", selected = "DEGs")
  })

  observeEvent(input$reset.selection1, {
    v$scData <- SetAllIdent(object = v$scData,  id = 'cluster.ident')
    #event_data("plotly_select") <- NULL
    v$plotlySelection <- NULL
  })

  observeEvent(input$PDFi, {
    if(!is.null(v$scData)){
      withProgress(message="Downloading plot PDF files...", value=0, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"UMAP_plot_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,"UMAP_plot_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        umapplot1 <- DimPlot(v$scData, reduction = "umap", label = T)
        umapplot2 <- DimPlot(v$scData, reduction = "umap", label = T, group.by = 'batch')
        umapplot3 <- DimPlot(v$scData, reduction = "umap", label = T, group.by = 'celltype')
        prePlot()
        pdf(filename2,
            width=as.numeric(input$pdf_w),
            height=as.numeric(input$pdf_h))
        print(umapplot1)
        print(umapplot2)
        print(umapplot3)
        dev.off()
      })
      withProgress(message="Downloading UMAP coordinates...", value=0.6, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"umap_", Sys.Date(), ".txt")
        j = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,"umap_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          j = j + 1;
        }
        write.csv(v$scData@reductions$umap@cell.embeddings, file = filename2)
      })
    }
  })

  ##---------------TSNE tabset-------------------
  output$perplex.option <- renderUI({
    if(is.null(v$isPCAdone)){
      return(NULL)
    }else{
      ##perplexity test
      n.cells <- isolate(nrow(v$scData@reductions$pca@cell.embeddings))
      max.perplex <- as.integer((n.cells - 1)/3)
      numericInput("perplexity",
                   label = "Perplexity",
                   value = if(max.perplex <30) max.perplex else 30,
                   min = 0,
                   max = max.perplex)
    }
  })

  observeEvent(input$doTsne, {
    withProgress(message = "Running tSNE...", value = 0.3, {
      #dims.use <- NULL
      #if(!is.null(v$scData@reductions$pca@jackstraw)){
      # score.df <- JackStraw_pval(v$scData)
      # dims.use <- signif_PCs(score.df)
      # }else{
      # dims.use <- 1:input$dim.used
      #}
      v$scData <- RunTSNE(v$scData, dims = 1:input$dim.used, perplexity = input$perplexity, assay = input$assays1)
      output$Tsne.done <- renderText(paste0("TSNE done!"))
      v$isTSNEdone <- TRUE
    })
  })

  output$Tsne_2d_plot_1 <- renderPlotly({
    if(is.null(v$scData) || is.null(v$isTSNEdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating TSNE 2D Plot...", value=0, {
        DimPlot(v$scData, reduction = "tsne", label = T)
      })
    }
  })

  output$Tsne_2d_plot_2 <- renderPlotly({
    if(is.null(v$scData) || is.null(v$isTSNEdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating TSNE 2D Plot...", value=0, {
        DimPlot(v$scData, reduction = "tsne", label = T, group.by = 'batch')
      })
    }
  })

  output$Tsne_2d_plot_3 <- renderPlotly({
    if(is.null(v$scData) || is.null(v$isTSNEdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating TSNE 2D Plot...", value=0, {
        DimPlot(v$scData, reduction = "tsne", label = T, group.by = 'celltype')
      })
    }
  })

  # output$Tsne_3d_plot <- renderPlotly({
  #     if(is.null(v$scData) || is.null(v$isTSNEdone)){
  #         plotly_empty()
  #    }else{
  #          withProgress(message="Generating TSNE 3D Plot...", value=0, {
  #             DimPlot(v$scData, reduction = "tsne", label = T)
  #          })
  #     }
  #  })

  output$selection.summary <- renderText({
    if(is.null(v$plotlySelection)){
      return(NULL)
    }else{
      t <- paste0(length(v$plotlySelection), " cells selected")
      t
    }
  })

  observeEvent(input$create.selection, {
    ## stash old identity
    if(is.null(v$scData@meta.data$cluster.ident)){
      v$scData <- StashIdent(object = v$scData, save.name = 'cluster.ident')
    }
    v$scData <- Idents(object = v$scData,
                       cells = v$plotlySelection,
                       value = as.character(input$selection.name)
    )
    updateTabsetPanel(session, "tabs", selected = "DEGs")
  })

  observeEvent(input$reset.selection, {
    v$scData <- SetAllIdent(object = v$scData,  id = 'cluster.ident')
    #event_data("plotly_select") <- NULL
    v$plotlySelection <- NULL
  })

  observeEvent(input$PDFj, {
    if(!is.null(v$scData)){
      withProgress(message="Downloading plot PDF files...", value=0, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"TSNE_plot_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,"TSNE_plot_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        tsneplot1 <- DimPlot(v$scData, reduction = "tsne", label = T)
        tsneplot2 <- DimPlot(v$scData, reduction = "tsne", label = T, group.by = 'batch')
        tsneplot3 <- DimPlot(v$scData, reduction = "tsne", label = T, group.by = 'celltype')
        prePlot()
        pdf(filename2,
            width=as.numeric(input$pdf_w),
            height=as.numeric(input$pdf_h))
        print(tsneplot1)
        print(tsneplot2)
        print(tsneplot3)
        dev.off()
      })
      withProgress(message="Downloading tSNE coordinates...", value=0.6, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"tsne_", Sys.Date(), ".txt")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,"tsne_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        write.csv(v$scData@reductions$tsne@cell.embeddings, file = filename2)
      })
    }
  })

  ##---------------DEGs tabset-------------------
  # output$clust1 <- renderUI({
  #  if(is.null(v$scData)){
  #    return(NULL)
  #  }else{
  #    celltypes <- levels(v$scData@active.ident)
  #     selectInput('c1', 'Choose cluster of interest:',
  #                choices = celltypes,
  #               selected = if("Selection" %in% celltypes) "Selection" else celltypes[1],
  #                selectize = FALSE,
  #                width = "100%")
  #  }
  # })
  # output$clust2 <- renderUI({
  #    if(is.null(v$scData)){
  #     return(NULL)
  #  }else{
  #    celltypes <- levels(v$scData@active.ident)
  #    selectInput('c2', 'Choose cluster to compare to:',
  #                choices = c("All", setdiff(celltypes, input$c1)),
  #                selected = "All",
  #               selectize = FALSE,
  #              width = "100%")
  #  }
  # })

  observeEvent(input$doDeg, {
    if(is.null(v$scData)){
      return(NULL)
    }else{
      withProgress(message="Finding DEGs...", value=0, {
        #      if(input$c2 == "All"){
        ips.markers <- FindAllMarkers(v$scData, only.pos = FALSE, min.pct = input$min_pct, logfc.threshold = input$logfc, assay = input$assays1)
        #     }else{
        #      ips.markers <- FindMarkers(v$scData, ident.1 = input$c1, ident.2 = input$c2, thresh.use = 2)
        #   }
        #  ips.markers$adj_p_val <- p.adjust(ips.markers$p_val, method = "BH")
        v$ips.markers <- ips.markers
      })
    }
  })

  output$deg.gene.select <- renderUI({
    if(is.null(v$ips.markers)){
      return(NULL)
    }else{
      selectInput("deg.gene", label = "Gene to visualise",
                  choices = rownames(v$ips.markers))
    }
  })

  output$Deg.plot <- renderPlotly({
    if(is.null(v$ips.markers)){
      return(NULL)
    }else{
      withProgress(message="Generating DEG Plot...", value=0, {
        VlnPlot(v$scData, input$deg.gene)
      })
    }
  })

  output$Deg1.plot <- renderPlotly({
    if(is.null(v$ips.markers)){
      return(NULL)
    }else{
      withProgress(message="Generating DEG Plot...", value=0, {
        FeaturePlot(v$scData, input$deg.gene)
      })
    }
  })

  output$Deg2.plot <- renderPlotly({
    if(is.null(v$ips.markers)){
      return(NULL)
    }else{
      withProgress(message="Generating DEG Plot...", value=0, {
        DotPlot(v$scData, features = input$deg.gene)
      })
    }
  })

  output$Deg3.plot <- renderPlotly({
    if(is.null(v$ips.markers)){
      return(NULL)
    }else{
      withProgress(message="Generating DEG Plot...", value=0, {
        v$ips.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
        DoHeatmap(v$scData, features = top10$gene, size = 4) + theme(axis.text.y = element_text(size = 5)) + NoLegend()
      })
    }
  })

  output$Deg.table <- DT::renderDataTable(
    #  if(is.null(v$scData)){
    #    return(NULL)
    #  }else{
    #    signif(v$ips.markers, 4)
    #  }
    v$ips.markers, options = list(scrollX = TRUE, scrollY = "400px"))

  observeEvent(input$PDFk, {
    if(!is.null(v$scData)){
      withProgress(message="Downloading plot PDF files...", value=0, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"DEG_plot_", input$c1, "vs", input$c2, "_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,"DEG_plot_", input$c1, "vs", input$c2, "_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        degVln <- VlnPlot(v$scData, input$deg.gene)
        degFeature <- FeaturePlot(v$scData, input$deg.gene)
        degDot <- DotPlot(v$scData, features = input$deg.gene)
        prePlot()
        pdf(filename2,
            width=as.numeric(input$pdf_w),
            height=as.numeric(input$pdf_h))
        print(degVln)
        print(degFeature)
        print(degDot)
        dev.off()
        write.csv(v$ips.markers, file = paste0(pdfDir, .Platform$file.sep,"DEG_table_", input$c1, "vs", input$c2, "_", Sys.Date(), ".csv"))
      })
    }
  })


  ##---------------Summary tab

  ##------Clean up when ending session----
  session$onSessionEnded(function(){
    prePlot()
  })
})


