source("library.R")$value

server <- function(input,output,session) 
{
  #Increase maximum allowable file size to 30 MB
  options(shiny.maxRequestSize=50*1024^2)
  
  # Adding code for rendering reports with code
  source("reports.R", local = TRUE)$value
  
  output$contents <-  DT::renderDataTable({
    data <- getNormalizedData()$raw
    if(!is.null(data)) data %>% DT::datatable(
      class = 'cell-border stripe',
      rownames = FALSE,
      filter   = 'top'
    )
  })
  
  output$downloadDEAFiles <- downloadHandler(
    filename = function() {
      paste("HOMRNA_DEA_results_files", "zip", sep=".")
    },
    content = function(fname) {
      fs <- c()
      data <- get.DEA.results()
      validate(
        need(!is.null(data), "Please, perform DEA analysis before exporting")
      )
      if (input$DEmethod == "deseq") {
        if (input$lfc) {
          resNames <- resultsNames(data)
          resNames <- resNames[-1]
          for (i in resNames){
            path <- paste0("ApeGLM",i,".csv")
            fs <- c(fs, path)
            res <- lfcShrink(data, type = "apeglm",
                             coef = i,
                             svalue=TRUE)
            tbl <-  cbind("Symbol"= rownames(res),as.data.frame(res))
            write_csv(tbl,path)
          }
        
      }
      
        else {
        for (i in resultsNames(data)){
          path <- paste0("DESeq",i,".csv")
          fs <- c(fs, path)
          df <- cbind("Symbol" = rownames(results(data,name = i)),as.data.frame(results(data,name = i)))
          write_csv(df, path)
          }
        }
      }

      else if (input$DEmethod =="edger") {
        path <- "edgerR.csv"
        fs <- c(fs, path)
        df <- cbind("Symbol" = rownames(data$table), data$table)
        write_csv(df, path)
      }
      zip(zipfile = fname, files = fs, flags = "-j")
    }, contentType = "application/zip"
  )  
  
  output$downloadNormalizedData <- downloadHandler(
    filename = function() {
      paste("HOMRNA_normalized_files", "zip", sep=".")
    },
    content = function(fname) {
      fs <- c()
      data <- getNormalizedData()
      for (i in names(data)){
        path <- paste0(i,".csv")
        fs <- c(fs, path)
        write.csv(data[[i]], path)
      }
      zip(zipfile=fname, files=fs, flags = "-j")
    },
    contentType = "application/zip"
  )  
  
  output$downloadPCA <- downloadHandler(
    filename = function() {
      paste("PCA_output_file", "csv", sep=".")
    },
    content = function(file) {
      data <- d
      write.csv(data, file)
    },
    contentType = "text/csv"
  )
  
  readData <- reactive({
    ret <- NULL
    inFile <- input$rawcounts
    if (!is.null(inFile))  {
      withProgress(message = 'Reading the data',
                   detail = "This may take a while", value = 0, {
                     ret <-  read_csv(inFile$datapath, col_types = readr::cols())
                     setProgress(1, detail = paste("Completed"))
                   }
      )
      if(!is.data.frame(ret)){
        withProgress(message = 'Reading the data',
                     detail = "This may take a while", value = 0, {
                       ret <-  read_csv2(inFile$datapath, col_types = readr::cols())
                       setProgress(1, detail = paste("Completed"))
                     }
        )
      }
    }  
    # Check if the input data has the required format
    ret <- checkDataInput(ret)
    ret
  })
  
  checkDataInput <- function(data){
    if(is.null(data)) return(NULL)
    if("status" %in% colnames(data)) data$status <- NULL
    
    tryCatch({
      res <- getEndGeneInfo(data)
    }, error = function(e){
      sendSweetAlert(
        session = session,
        title =  "Input error", 
        text =   paste0("Test1Data uploaded does not have the expected format.\n", 
                        "\nWe were unable to identify the gene column or map it to hg38 or mm10.", 
                        "\nThe expected input is a table with genes in the first column and genes raw counts on the other ones."),
        type = "error"
      )
      return(NULL)
    })
#    if(!"Chr" %in% colnames(res$data)){
#      sendSweetAlert(
#        session = session,
#        title =  "Input error", 
#        text =   paste0("Test2Data uploaded does not have the expected format.\n", 
#                        "\nWe were unable to identify the gene column or map it to hg38 or mm10.", 
#                        "\nThe expected input is a table with genes in the first column and genes raw counts on the other ones."),
#        type = "error"
#      )
#      return(NULL)
#    }
    
    data <- tryCatch({
      colSums(res$data[,(res$ngene + 1):ncol(res$data)])
      return(data)
    }, error = function(e){
      
      sendSweetAlert(
        session = session,
        title =  "Input error", 
        text =   paste0("Test3Data uploaded does not have the expected format.\n", 
                        "\nWe were unable to identify the gene column or map it to hg38 or mm10.", 
                        "\nThe expected input is a table with genes in the first column and genes raw counts on the other ones."),
        type = "error"
      )
      return(NULL)
    })
    return(data)
  }
  
  readMetaData <- reactive({
    ret <- NULL
    inFile <- input$metadata
    if (!is.null(inFile))  {
      withProgress(message = 'Reading the data',
                   detail = "This may take a while", value = 0, {
                     ret <-  read_csv(inFile$datapath, col_types = readr::cols())
                     setProgress(1, detail = paste("Completed"))
                   }
      )
    }
    ret
  })
  
  vennCounts <- function(group) {
    log2Cutoff <- input$vennLog2FoldChange
    pvalCutoff <- input$vennPadj
    data <- get.DEG.results()
    DEG1 <- NULL
    DEG2 <- NULL
    DEG3 <- NULL
    DEG4 <- NULL
    if(is.null(data)) return(NULL)
    if(!is.null(data$DEG1)) {
      DEG1 <- subset(data$DEG1, svalue<pvalCutoff,
                     select=c(Symbol,log2FoldChange))
      DEG1 <- subset(DEG1, log2FoldChange > log2Cutoff | log2FoldChange < -log2Cutoff,
                     select=Symbol)
      DEG1$Group <- input$DEfile1Name
    }
    if(!is.null(data$DEG2)) {
      DEG2 <- subset(data$DEG2, svalue<pvalCutoff,
                     select=c(Symbol,log2FoldChange))
      DEG2 <- subset(DEG2, log2FoldChange > log2Cutoff | log2FoldChange < -log2Cutoff,
                     select=Symbol)
      DEG2$Group <- input$DEfile2Name
    }
    if(!is.null(data$DEG3)) {
      DEG3 <- subset(data$DEG3, svalue<pvalCutoff,
                     select=c(Symbol,log2FoldChange))
      DEG3 <- subset(DEG3, log2FoldChange > log2Cutoff | log2FoldChange < -log2Cutoff,
                     select=Symbol)
      DEG3$Group <- input$DEfile3Name
    }
    if(!is.null(data$DEG4)) {
      DEG4 <- subset(data$DEG4, svalue<pvalCutoff,
                     select=c(Symbol,log2FoldChange))
      DEG4 <- subset(DEG4, log2FoldChange > log2Cutoff | log2FoldChange < -log2Cutoff,
                     select=Symbol)
      DEG4$Group <- input$DEfile4Name
    }
    allData <- rbind(DEG1, DEG2, DEG3,DEG4)
    allData <- drop_na(allData)
    ppl <- allData
    print("ppl worked")
    print(ppl)
    for (i in 1:length(group)) {
      ppl <- subset(ppl, Group == group[i])
    }
    nrow(ppl)
    print("nrow worked")
    print(nrow(ppl))
  }
  
  plotVenn <- function(groups, ...) {
#    grid.newpage()
    if (length(groups) == 1) {
      out <- draw.single.venn(vennCounts(groups), ...)
    }
    if (length(groups) == 2) {
      out <- draw.pairwise.venn(vennCounts(groups[1]), vennCounts(groups[2]), vennCounts(groups[1:2]), ...)
    }
    if (length(groups) == 3) {
      out <- draw.triple.venn(vennCounts(groups[1]), vennCounts(groups[2]), vennCounts(groups[3]), vennCounts(groups[1:2]), 
                              vennCounts(groups[2:3]), vennCounts(groups[c(1, 3)]), vennCounts(groups), ...)
    }
    if (length(groups) == 4) {
      out <- draw.quad.venn(vennCounts(groups[1]), vennCounts(groups[2]), vennCounts(groups[3]), vennCounts(groups[4]), 
                            vennCounts(groups[1:2]), vennCounts(groups[c(1, 3)]), vennCounts(groups[c(1, 4)]), vennCounts(groups[2:3]), 
                            vennCounts(groups[c(2, 4)]), vennCounts(groups[3:4]), vennCounts(groups[1:3]), vennCounts(groups[c(1, 2, 
                                                                                       4)]), vennCounts(groups[c(1, 3, 4)]), vennCounts(groups[2:4]), vennCounts(groups), ...)
    }
    return(out)
  }
  
  ### reactive fct that calcs the transforms and saves them so it doesnt take too long each time
  getNormalizedData <- reactive({
    if (!is.null(readData())) {
      all_cell_lines <- readData()
      
      # Add gene metadata information
      withProgress(message = 'Adding gene metadata',
                   detail = "This may take a while", value = 0, {
                     # We will check if metadata was added
                     tryCatch({
                       res <- getEndGeneInfo(all_cell_lines)
                       all_cell_lines <- res$data
                       ngene <- res$ngene
                     } , error = function(e){
                       sendSweetAlert(
                         session = session,
                         title =  "Error reading the data",
                         text =   paste0("Please check"),
                         type = "error"
                       )                     })
                   }
      )
      
      tbl.tab1 <- all_cell_lines[rowSums(all_cell_lines[,(ngene + 1):ncol(all_cell_lines)]) > 1,] ##filtering step, actually change the object
      data <- as.matrix(tbl.tab1[,(ngene + 1):ncol(tbl.tab1)])
      metadata <- tbl.tab1[,1:(ngene)]
      
      withProgress(message = 'Normalizing data',
                   detail = "This may take a while", value = 0, {
                     # normalization: rlog takes a lot of time (hours for a big matrix)
                     raw      <- cbind(metadata, data) # might might have to take out blind option???
                     setProgress(0.1, detail = paste("Starting VST"))
                     vst      <- cbind(metadata, tryCatch(vst(data),error = function(e){varianceStabilizingTransformation(data)}))
                     setProgress(0.2, detail = paste("VST completed, starting rownorm"))
                     rownorm  <- cbind(metadata, rownorm(data))
                     setProgress(0.5, detail = paste("rownorm completed, starting CPM"))
                     cpm      <- cbind(metadata, cpm(data, log = TRUE)) ##adding log=TRUE option
                     setProgress(0.7, detail = paste("CPM completed, starting rlog"))
                     if(ncol(data) < 30){
                       rlog     <- cbind(metadata, rlog(data))
                       ret      <- list(vst,rownorm,raw,cpm,rlog)
                       names(ret) <- c("vst","rownorm","raw","cpm","rlog")
                     } else {
                       ret      <- list(vst,rownorm,raw,cpm)
                       names(ret) <- c("vst","rownorm","raw","cpm")
                     }
                     setProgress(1, detail = paste("Completed"))
                   }
      )
      #save(ret,file = "genavi.rda")
    } else { ##### not sure if need to remove? leave in -> no change, comment out -> breaks app
      ret <- get(load("HomRNA.rda"))
    }
    if(length(names(ret)) == 5){
      transforms <-  c("raw counts",  
                       "row normalized",  
                       "logCPM - log Counts per Million",  
                       "vst - Variance Stabilizing Transformation",  
                       "rlog - regularized logarithm")
    } else {
      transforms <-  c("raw counts",  
                       "row normalized",  
                       "logCPM - log Counts per Million",  
                       "vst - Variance Stabilizing Transformation")
    }
    updateSelectizeInput(session, 'select_tab1', 
                         selected = "raw counts",
                         choices =  transforms,
                         server = TRUE)
    return(ret)
  })
  getTab1 <- reactive({
    data <- getNormalizedData()
    select <- input$select_tab1
    tbl.tab1 <- NULL
    if(select == "raw counts") tbl.tab1 <- data$raw #table.counts #DT::datatable(table.counts)
    if(select == "rlog - regularized logarithm")  tbl.tab1 <- data$rlog ##include warning in vignette, dependent on number of columns
    if(select == "vst - Variance Stabilizing Transformation")  tbl.tab1 <- data$vst
    if(select == "row normalized")  tbl.tab1 <- data$rownorm
    if(select == "logCPM - log Counts per Million")  tbl.tab1 <- data$cpm
    
    sampFile <- input$input_sample_list_tab1
    if (!is.null(sampFile)) {
      selections <- read_lines(sampFile$datapath)
#      print(sampList)
#      selections <- unique(c(selected_rows,which(matrix_clus[,1] %in% sampList)))
#      print(selections)
    }
    else {selections <- input$selectedSamples
    print(selections)
    }
    
    if(!is.null(selections)) {tbl.tab1 <- tbl.tab1 %>% select(c("Symbol",selections))
    }
    else {tbl.tab1 <- tbl.tab1
    }
    tbl.tab1
  })
  
  
  gene.selection <- reactiveVal(0)      
  observeEvent(input$select_most_variable, {
    gene.selection("Most Var")             # rv$value <- newValue
  })
  observeEvent(input$unselect_all, {
    gene.selection("Unselect")             # rv$value <- newValue
  })
  
  observeEvent(input$input_gene_list_but, {
    aux <- gene.selection()
    if(is.character(aux)) aux <- 0
    gene.selection( aux + 1)             # rv$value <- newValue
  })
  
  sample.selection <-reactiveVal(0)
  observeEvent(input$selectedSamples, {
    aux <- sample.selection()
    if(is.character(aux)) aux <- 0
    gene.selection( aux + 1)             # rv$value <- newValue
  })
  
  
  # Function to add which genes were selected by the user
  getTab1.selected <- reactive({
    
    # get original table
    tab <- getTab1()
    if(gene.selection() == "Most Var"){
      m <- as.matrix(tab[,grep("EnsemblID|Symbol|Genename|Geneid|Chr|Start|End|Strand|Length", colnames(tab),ignore.case = T,invert = T)])
      selected_rows <- sort(rowVars(m,na.rm = T), decreasing = T,index.return=TRUE)$ix[1:1000]
      status <- factor("Unselected",levels = c("Unselected","Selected"))
      tab <- cbind(status,tab)
      tab$status[selected_rows] <- 'Selected'
    } else if (gene.selection() == "Unselect") {
      status <- factor("Unselected",levels = c("Unselected","Selected"))
      tab <- cbind(status,tab)
      selected_rows <- NULL
    } else {
      # which ones are already selected
      selected_rows <- input$tbl.tab1_rows_selected
      
      # which column has our gene symbol ?
      idx <- grep("symbol|genename",colnames(tab),ignore.case = T)
      
      # Read genes to filter from a file 
      inFile <- input$input_gene_list_tab1
      if (!is.null(inFile)) {
        geneList <- read_lines(inFile$datapath)
        selected_rows <- unique(c(selected_rows,which(tab[,idx] %in% geneList)))
      }    
      
      # Read genes to filter from a textarea 
      geneList <- isolate({input$input_gene_list_area})
      if(!is.null(geneList)){
        geneList <- parse.textarea.input(geneList)
        selected_rows <- unique(c(selected_rows,which(tab[,idx] %in% geneList)))
      }
      
      # Update the rows selected
      status <- factor("Unselected",levels = c("Unselected","Selected"))
      tab <- cbind(status,tab)
      tab$status[selected_rows] <- 'Selected'
    }
    return(list("tab" = tab,"selected_rows" = selected_rows))
    
  })
  
  output$tbl.tab1 <-  DT::renderDataTable({
    tbl.tab1 <- getTab1.selected()$tab
    selected_rows <- getTab1.selected()$selected_rows
    if(is.null(tbl.tab1)) return(NULL)
    ######### sorting by mean and sd ##################### ....fucks up the select sorting thing...
    #if(input$select_sort_tab1 == "-no selection-") {return(tbl.tab1)}
    #if(input$select_sort_tab1 == "mean")
    #{
    # tbl.tab1 <- tbl.tab1[order(apply(tbl.tab1[,8:dim(tbl.tab1)[2]],1,mean), decreasing = TRUE),]
    #}
    #if(input$select_sort_tab1 == "standard deviation")
    #{
    # tbl.tab1 <- tbl.tab1[order(apply(tbl.tab1[,8:dim(tbl.tab1)[2]],1,sd), decreasing = TRUE),]
    #}
    
    ######## this section sorts the table so that selected rows are first #######
    ####### ordering rows like this makes the selection wonky in the figures, other rows than what you select are being displayed
    ####### try taking this out to see how tables are rendered...or not rendered???
    tbl.tab1 %>% createTable(selected_rows,tableType = isolate({input$select_tab1}))
    
    ## try adding the genes list to match() here, see if it breaks the app
    ##gene_list_tab1 <- input$input_gene_list_tab1 ##create object here see if it breaks the app.....yup breaks the app
    
  }) ##works to get selected rows on top but fucks up if select more than one at a time.....worry about it later
  
  sel <- reactive({!is.null(input$tbl.tab1_rows_selected)})  
  
  observeEvent(sel(), {
    closeAlert(session, "geneAlert2")
    
    if(!sel()){
      createAlert(session, "genemessage", "geneAlert2", title = "Missing data", style =  "danger",
                  content = paste0("Please select genes in Gene Expression tab"),
                  append = FALSE)
      shinyjs::hide(id = "expression_heatmap", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
      shinyjs::hide(id = "expression_plots", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    } else if(length(input$tbl.tab1_rows_selected) == 1){
      shinyjs::show(id = "expression_plots", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
      shinyjs::hide(id = "expression_heatmap", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    } else {
      shinyjs::show(id = "expression_heatmap", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
      shinyjs::hide(id = "expression_plots", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    }
  })
  
  output$barplot <- renderPlotly({
    
    if(is.null(input$tbl.tab1_rows_selected)) {
      createAlert(session, "genemessage", "geneAlert2", title = "Missing data", style =  "danger",
                  content = paste0("Please select genes in Gene Expression tab"),
                  append = FALSE)
      return(NULL) 
    }
    if(length(input$tbl.tab1_rows_selected) == 0) {
      createAlert(session, "genemessage", "geneAlert2", title = "Missing data", style =  "danger",
                  content = paste0("Please select genes in Gene Expression tab"),
                  append = FALSE)
      return(NULL) 
    }
    ##may need to put this in heatmap section and in tab2
    if(length(input$tbl.tab1_rows_selected) > 1) return(NULL)
    tbl.tab1 <- getTab1()
    
    # Columns 1 to 7: Genename  Geneid Chr   Start   End Strand Length 
    res <- getEndGeneInfo(tbl.tab1)
    ngene <- res$ngene
    
    tbl.tab1 <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected)
    p <- as.data.frame(t(tbl.tab1[,(res$ngene+1):ncol(tbl.tab1)]))
    colnames(p) <- "value"
    p$cell_line <- rownames(p)
    order <- rownames(p)
    
    
    select <- isolate(input$select_tab1)
    barplot <- ggplot(p, aes(x=cell_line, y=value)) + 
      geom_bar(stat = "identity") +  
      theme_bw() + 
      labs(x =  "", y = select) +
      scale_x_discrete(limits = order) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    ggplotly(barplot)
  })
  
  #---------------------------------
  # Heatmap plot tab
  #---------------------------------
  # calculate the variance for each gene
  output$pca_plot <- renderPlotly({
    nb.cols <- 21
    mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
    tbl.tab1 <- getTab1()

    
    # Columns 1 to 7: Genename  Geneid Chr   Start   End Strand Length 
    res <- getEndGeneInfo(tbl.tab1)
    ngene <- res$ngene
    
    m <- tbl.tab1 %>% dplyr::select((res$ngene + 1):ncol(tbl.tab1)) %>% as.matrix 
    
    df <- tbl.tab1 #%>% dplyr::select((res$ngene + 1):ncol(tbl.tab1)) 
    df <- data.frame(tbl.tab1)


    select <- 1:nrow(m)
    if(input$select_pca_type == "Top 1000 variable genes"){
      ntop <- 1000
      rv <- rowVars(m)
      
      # select the ntop genes by variance
      select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    } else   if(input$select_pca_type == "Selected genes"){
      select <- input$tbl.tab1_rows_selected
      if(length(select) < 3) {
        sendSweetAlert(
          session = session,
          title =  "Missing data",
          text =   paste0("Please select at least 3 genes in Gene Expression tab"),
          type = "info"
        )
        return(NULL)
      }
    }
    
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(m[select,]))
#    pca <- PcaGrid(t(m[select,]))
#    print(pca)
    df <- df[select,]
    rownames(df) <- df$Symbol
    df <- df[ , -which(names(df) %in% c("Symbol"))]

    pca2 <- prcomp(t(df))

    test <- mvnorm.skew.test(m[select,])

    print(test)
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    if(!is.null(input$pcacolor) & stringr::str_length(input$pcacolor) > 1) {
      sample.id <- which(readMetaData()[1,] %in% colnames(m))[1]
      color.idx <- which(input$pcacolor == colnames(readMetaData()))
      metadata <- readMetaData()[,c(1,color.idx)]
      d <<- data.frame(PC1 = pca$x[,1], 
                      PC2 = pca$x[,2], 
                      PC3 = pca$x[,3], 
                      name = colnames(m),
                      color = metadata[match(colnames(m),metadata[,1,drop = T]),2,drop = T])
      #This uses the ellipse function to generate 95% confidence ellipse on plots
      colors <- metadata[match(colnames(m),metadata[,1,drop = T]),2,drop = T]

    } else {
      d <<- data.frame(PC1 = pca$x[,1], 
                      PC2 = pca$x[,2], 
                      PC3 = pca$x[,3], 
                      name = colnames(m),
                      color = rep("blue",length(colnames(m))))
    }
#    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    if(input$pca_dimensions == "2D" & input$confEllipse == "confEll") {
      ind.p <- fviz_pca_ind(pca,
                        geom.ind="point",
                        pointshape=20,
                        pointsize=1.5,# show points only (nbut not "text")
                        col.ind = metadata[match(colnames(m),metadata[,1,drop = T]),2,drop = T], # color by groups
                        addEllipses = TRUE,
                        ellipse.type="norm",
                        ellipse.level=0.95,
                        legend.title = "Groups",
                        palette="mycolors",
                        mean.point=FALSE,
                        xlab = paste("PC1: ", round(percentVar[1] * 100, digits = 2), "% variance", sep=""),
                        ylab = paste("PC2: ", round(percentVar[2] * 100, digits = 2),"% variance", sep=""))
    } 
    
    else if(input$pca_dimensions == "2D" & input$confEllipse == "noConf" & input$geneLoadings == "noGeneLoad") {
      ind.p <- fviz_pca_ind(pca,
                            geom.ind="point",
                            pointshape=20,
                            pointsize=1.5,# show points only (nbut not "text")
                            col.ind = metadata[match(colnames(m),metadata[,1,drop = T]),2,drop = T], # color by groups
                            legend.title = "Groups",
                            palette="mycolors",
                            mean.point=FALSE,
                            xlab = paste("PC1: ", round(percentVar[1] * 100, digits = 2), "% variance", sep=""),
                            ylab = paste("PC2: ", round(percentVar[2] * 100, digits = 2),"% variance", sep=""))
    }
    

    else if(input$pca_dimensions == "2D" & input$confEllipse == "noConf" & input$geneLoadings == "GeneLoad") {
      ind.p <- fviz_pca_biplot(pca2,
                            geom.ind="point",
                            pointshape=20,
                            pointsize=1.5,# show points only (nbut not "text")
                            col.ind = metadata[match(colnames(m),metadata[,1,drop = T]),2,drop = T], # color by groups
                            legend.title = "Groups",
                            palette="mycolors",
                            mean.point=FALSE,
                            xlab = paste("PC1: ", round(percentVar[1] * 100, digits = 2), "% variance", sep=""),
                            ylab = paste("PC2: ", round(percentVar[2] * 100, digits = 2),"% variance", sep=""),
                            select.var=list(contrib=input$topGenes))
      
      pca3 <- FactoMineR::PCA(t(df))
      contrib_PC12 = pca3$var$contrib[,1:2]  # get contributions of each variable to PC1 and PC2
      eig_PC12 = pca3$eig[,1][1:2] # get the eigenvalues of both components
      contrib_total = apply(contrib_PC12, 1, function(x) {sum(x*eig_PC12)/sum(eig_PC12)}) # calculate the total contribution
      N = input$topGenes
      names(sort(contrib_total, decreasing = TRUE)[1:N]) -> topN_PC12
      tbl.tab3 <- topN_PC12
    }
    else {
      p <-
        plot_ly(
          d,
          x = ~ PC1 ,
          y = ~ PC2,
          z = ~ PC3,
          color = ~ color,
          text = ~ paste(name),
          type = "scatter3d",
          marker = list(size = 14),
          width = 800,
          height = 600
        ) %>%
        add_markers()
      ind.p <- layout(p, 
                  scene = list(
                    title = "Principal Component Analysis", 
                    xaxis = list(title = paste0("PC1: ", round(percentVar[1] * 100, digits = 2),"% variance")), 
                    yaxis = list(title = paste0("PC2: ", round(percentVar[2] * 100, digits = 2),"% variance")),
                    zaxis = list(title = paste0("PC3: ", round(percentVar[3] * 100, digits = 2),"% variance")) 
                  )
      )
      if(is.null(input$pcacolor) || stringr::str_length(input$pcacolor) == 0) { 
        ind.p <- layout(p,showlegend = FALSE)
      }
    }
    
    ind.p
  })
  
  output$heatmap_expr <- renderIheatmap({


    # we can't do a heatmap with only one gene
    if(length(input$tbl.tab1_rows_selected) < 2) return(NULL)
    
    tbl.tab1 <- getTab1()
    
    
    if("Symbol" %in% colnames(tbl.tab1)){
      geneNames <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected) %>% pull("Symbol")
    } else if("Genename" %in% colnames(tbl.tab1))  {
      geneNames <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected) %>% pull("Genename")
    } else {
      geneNames <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected) %>% pull(1)
    }
    
    # Columns 1 to 7: Genename  Geneid Chr   Start   End Strand Length 
    res <- getEndGeneInfo(tbl.tab1)
    ngene <- res$ngene
    color.idx <- which(input$Expheatcolor == colnames(readMetaData()))
    metadata <- readMetaData()[,c(1,color.idx)]

    matrix_expr <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected) %>% dplyr::select((res$ngene+1):ncol(tbl.tab1)) %>% as.matrix
    metadata <- readMetaData()[,c(1:8)] #Uncomment this to use all columns 2-8 of metadat for heatmap
    forAnnot <- metadata[match(colnames(matrix_expr), metadata$Sample),] #Comment this out to use the above line
#    forAnnot <- subset(metadata, select=c(input$Expheatcolor)) #Comment this out to use the alternative forAnnot above
    
    font.size <- ifelse(ncol(matrix_expr) > 30, 6, 12)
    font.size.genes <- ifelse(nrow(matrix_expr) > 30, 6, 12)
    
    heatmap_expr <- main_heatmap(matrix_expr,  name = isolate(input$select_tab1)) %>%
      add_col_annotation(forAnnot) %>%
      add_col_labels(ticktext = colnames(matrix_expr),font = list(size = font.size),size = 0.2) %>%
      add_row_labels(ticktext = geneNames,font = list(size = font.size),size = 0.2) %>% 
      add_col_dendro(hclust(dist(t(as.matrix(matrix_expr)))), reorder = TRUE) 
    
    if(nrow(matrix_expr) > 1) 
    {
      heatmap_expr <- heatmap_expr %>% add_row_dendro(hclust(dist((as.matrix(matrix_expr)))), reorder = TRUE, side = "right") 
    } 
    heatmap_expr
  })
  
  
  output$heatmap_clus <- renderIheatmap({
    closeAlert(session, "geneAlert")
    tbl.tab2 <- getTab1()
    
    # Columns 1 to 7: Genename  Geneid Chr   Start   End Strand Length 
    res <- getEndGeneInfo(tbl.tab2)
    ngene <- res$ngene
    
    matrix_clus <- tbl.tab2[,c(1,(res$ngene+1):ncol(tbl.tab2))] 
    
    #replace above command with this based on select input
    if(input$select_clus == "-no selection-") return(NULL) ##commenting it out still has filtered hm show automatically
    ##BT549 disapears from list of cell lines???
    ##how to make this heatmap show by default/automatically
    ##doesn't work with raw counts
    if(input$select_clus == "All genes") {
      data <- matrix_clus[,-1]
    } else { 
      # get selected genes
      selected_rows <- input$tbl.tab1_rows_selected
      if(length(selected_rows) < 1) {
        sendSweetAlert(
          session = session,
          title =  "Missing data",
          text =   paste0("Please select genes in Gene Expression tab"),
          type = "error"
        )
        return(NULL)
      }
      inFile <- input$input_gene_list_tab1
      if (!is.null(inFile)) {
        geneList <- read_lines(inFile$datapath)
        selected_rows <- unique(c(selected_rows,which(matrix_clus[,1] %in% geneList)))
      }
      data <- matrix_clus[selected_rows,-1]
    }
    
    # The deafult is to cluster by sample, but we will cluster by gene if selected
    if(input$select_clus_type == "Genes") {
      selected_rows <- input$tbl.tab1_rows_selected
      if(length(selected_rows) < 1) {
        sendSweetAlert(
          session = session,
          title =  "Missing data",
          text =   paste0("Please select genes in Gene Expression tab"),
          type = "error"
        )
        return(NULL)
      }
      inFile <- input$input_gene_list_tab1
      if (!is.null(inFile)) {
        geneList <- read_lines(inFile$datapath)
        selected_rows <- unique(c(selected_rows,which(matrix_clus[,1] %in% geneList)))
      }
      data <- matrix_clus[selected_rows,-1]
      rownames(data) <- matrix_clus[selected_rows,1]
      data <- t(data)
    }
    
    font.size <- ifelse(ncol(data) > 30, 6, 12)
    heatmap_clus <-  tryCatch({
      main_heatmap(as.matrix(cor(data, method = "pearson")), name = "Correlation", colors = custom_pal_blues) %>% 
        add_col_labels(ticktext = colnames(data),font = list(size = font.size),size = 0.2) %>%
        add_row_labels(ticktext = colnames(data),font = list(size = font.size),size = 0.2) %>% 
        add_col_dendro(hclust(as.dist(1 - cor(data, method = "pearson"))), reorder = TRUE) %>% 
        add_row_dendro(hclust(as.dist(1 - cor(data, method = "pearson"))), reorder = TRUE, side = "right") 
    }, warning = function(w){
      sendSweetAlert(
        session = session,
        title =  "Sorry, we had an error...",
        text =  paste0("Clustering is not possible.",
                       "We tried to cluster ", input$select_clus_type ," using ", input$select_clus,"."),
        type = "error"
      )
      return(NULL)
    }, error = function(e){
      
      sendSweetAlert(
        session = session,
        title =  "Sorry, we had an error...",
        text =  "Clustering is not possible",
        type = "error"
      )
      return(NULL)
    })
    if(is.null(heatmap_clus)) return(NULL)
    heatmap_clus
  })
  
  #------------------------------------------
  # DEA - differential expression analysis
  #------------------------------------------
  output$downloadData <- downloadHandler(
    filename = function() {
      "Chlam_metaData.csv"
    },
    content = function(file) {
      metadata <- readr::read_csv("Chlam_metaData.csv")
      write_csv(metadata, file)
    }
  )
  
  observeEvent(input$select_clus_type,{
    if(input$select_clus_type == "Sample") choices <- c("All genes", "Selected genes")
    if(input$select_clus_type == "Genes") choices <- c("Selected genes")
    updateSelectizeInput(session, 'select_clus', choices = choices, server = TRUE)
  })  
  observe({
    metadata <- readMetaData() ### try reactive?
    if(!is.null(metadata)){
      updateSelectizeInput(session, 'condition', choices =  colnames(metadata)[-1], server = TRUE)
    }
    if(!is.null(metadata)){
      updateSelectizeInput(session, 'covariates', choices =  c(" ",colnames(metadata)[-1]), server = TRUE)
    }
    if(!is.null(metadata)){
      updateSelectizeInput(session, 'pcacolor', choices =  colnames(metadata)[-1], server = TRUE)
      shinyjs::show("pcacolor")
      updateSelectizeInput(session, 'Expheatcolor', choices =  colnames(metadata)[-1], server = TRUE)
      shinyjs::show("Expheatcolor")
      data <- getNormalizedData()
      updateSelectizeInput(session, 'selectedSamples', choices =  colnames(data$raw)[-1], server = TRUE) 
      shinyjs::show("selectedSamples")
    } else {
      shinyjs::hide("pcacolor")
      shinyjs::hide("Expheatcolor")
      data <- getNormalizedData()
      updateSelectizeInput(session, 'selectedSamples', choices =  colnames(data$raw)[-1], server = TRUE)#This is included here in case we want to hide this based on the presence of metadata
      shinyjs::show("selectedSamples")
    }
  })
  
  get.DEA.results <- reactive({
    input$dea
    metadata <- readMetaData()
    if(is.null(metadata)) {
      sendSweetAlert(
        session = session,
        title =  "Missing metadata",
        text = paste0("Please upload metadata file"),
        type = "error"
      )
      return(NULL)
    }
    if(!is.null(readData())) all_cell_lines <- readData()
    
    res <- getEndGeneInfo(all_cell_lines)
    all_cell_lines <- res$data
    ngene <- res$ngene
    
    if("Symbol" %in% colnames(all_cell_lines)){
      genes <- all_cell_lines %>% pull("Symbol")
    } else if("Genename" %in% colnames(all_cell_lines))  {
      genes <- all_cell_lines %>% pull("Genename")
    } else {
      genes <- all_cell_lines %>% pull(1)
    }
    
    cts <- as.matrix(all_cell_lines[,(ngene + 1):ncol(all_cell_lines)])
    rownames(cts) <-  genes
    
    # Read aux values required for analysis (condition, covariates and reference value)
    cond <- isolate(input$condition)
    cov <- isolate(input$covariates)
    ref <-  isolate(input$reference)
    DEmeth <- isolate(input$DEmethod)
    
    if(is.null(cond))   {
      sendSweetAlert(
        session = session,
        title =  "Missing condition",
        text = paste0("Please select condition file"),
        type = "error"
      )
      return(NULL)
    } 
    form <- getFormula()
    if(is.null(form)){
      sendSweetAlert(
        session = session,
        title =  "Missing formula",
        text = paste0("Please select condition column"),
        type = "error"
      )
      return(NULL)
    }
    if(str_length(ref) == 0)   {
      
      sendSweetAlert(
        session = session,
        title =  "Missing reference level",
        text = paste0("Please select reference level"),
        type = "error"
      )
      return(NULL)
    } 
    if(nrow(metadata) != ncol(cts))   {
      sendSweetAlert(
        session = session,
        title =  "Metadata error",
        text = paste0("Metadata and data does not have same samples"),
        type = "error"
      )
      return(NULL)
    } 
    
    withProgress(message = 'DEA Analysis',
                 detail = "Creating input file", value = 0, {
                   
                   if(!all(metadata %>% pull(1) %in% colnames(cts)))   {
                     sendSweetAlert(
                       session = session,
                       title =  "Metadata error",
                       text = paste0("First column of the metadata file must have the mapping to the samples with the exact names"),
                       type = "error"
                     )
                     return(NULL)
                   } 
                   metadata <- metadata[match(colnames(cts), metadata %>% pull(1)),]
                   if(!all(metadata %>% pull(1) == colnames(cts)))   {
                     
                     sendSweetAlert(
                       session = session,
                       title =  "Metadata error",
                       text = paste0("First column of the metadata file must have the mapping to the samples with the exact names"),
                       type = "error"
                     )
                     return(NULL)
                   } 
                   dds <-  tryCatch({
                     sampDEAFile <- input$input_sample_list_DEA
                     if (!is.null(sampDEAFile)) {
                       keep.samples <- read_lines(sampDEAFile$datapath)
                       newMetaData <- metadata %>% filter(Samples %in% keep.samples)
                       newCTS <- cts[ , keep.samples]
                     } else {
                       newMetaData <- metadata
                       newCTS <- cts
                     }
                       
                     keep.samples <- !is.na(newMetaData[,input$condition,drop = T])
                     print(length(keep.samples))
                     
                     if(any(is.na(newMetaData[,input$condition,drop = T]))){
                       
                       sendSweetAlert(
                         session = session,
                         title =  paste0(sum(!keep.samples), " samples with have NA annotations"), 
                         text = "To perform the DEA samples cannot be labled as NA we will remove it",
                         type = "error"
                       )
                     }
                     
                     if (DEmeth=="deseq") {
                       dds <- DESeqDataSetFromMatrix(countData = newCTS[,keep.samples],
                                                     colData = newMetaData[keep.samples,],
                                                     design = form)
                       setProgress(0.2, detail = paste("Performing DEA DESeq2"))
                       keep <- rowSums(counts(dds)) >= 10
                       dds[[cond]] <- relevel(dds[[cond]], ref = ref)
                       dds <- dds[keep,]
                       dds <- DESeq(dds)
                       return(dds)
                     }
                     
                     else if (DEmeth=="edger") {
                       keep.samples <- !is.na(metadata[,input$condition,drop = T])
                       setProgress(0.2, detail = paste("Performing DEA edgeR-QL"))
                       keepMeta <- metadata[keep.samples,]
                       group <- factor(keepMeta[[cond]])
                       edgerY <- DGEList(counts=cts, group=group)
                       keep.samps <- filterByExpr(edgerY, group=group)
                       edgerY <- edgerY[keep.samps,, keep.lib.sizes=FALSE]
                       edgerY <- calcNormFactors(edgerY)
                       design <- model.matrix(~0+group)
                       colnames(design) <- levels(group)
                       conts <- edgerConts(colnames(design),ref)
                       edgerY <- estimateDisp(edgerY, design, robust=TRUE)
                       dds <- glmQLFit(edgerY, design, robust=TRUE)
                       newConts <-c()
                       for (i in conts) {
                         newConts <- append(newConts,i)
                       }
                       my.contrasts <- makeContrasts(contrasts=newConts,levels=design)
                       dds <- glmQLFTest(dds, contrast=my.contrasts)
                       return(dds)
                     }
                   }, error = function(e){
                     
                     sendSweetAlert(
                       session = session,
                       title = "Error in DEA",
                       text = paste0(e),
                       type = "error"
                     )
                     
                     return(NULL)
                   })
                 }
    )
    
    return(dds)
  })
  
  get.DEG.results <- reactive({
    ret1 <- NULL
    ret2 <- NULL
    ret3 <- NULL
    ret4 <- NULL
    ret5 <- NULL
    a <- 0
    b <- 0
    c <- 0
    d <- 0
    e <- 0
    inFile1 <- input$DEfile1
    inFile2 <- input$DEfile2
    inFile3 <- input$DEfile3
    inFile4 <- input$DEfile4
    inFile5 <- input$DEfile5
    if (!is.null(inFile1))  {
      withProgress(message = 'Reading the data',
                   detail = "This may take a while", value = 0, {
                     ret1 <-  read_csv(inFile1$datapath, col_types = readr::cols())
                     a <- 1
                     setProgress(1, detail = paste("Completed"))
                   }
      )
    }
    if (!is.null(inFile2))  {
      withProgress(message = 'Reading the data',
                   detail = "This may take a while", value = 0, {
                     ret2 <-  read_csv(inFile2$datapath, col_types = readr::cols())
                     b <- 1
                     setProgress(1, detail = paste("Completed"))
                   }
      )
    }
    if (!is.null(inFile3))  {
      withProgress(message = 'Reading the data',
                   detail = "This may take a while", value = 0, {
                     ret3 <-  read_csv(inFile3$datapath, col_types = readr::cols())
                     c <- 1
                     setProgress(1, detail = paste("Completed"))
                   }
      )
    }
    if (!is.null(inFile4))  {
      withProgress(message = 'Reading the data',
                   detail = "This may take a while", value = 0, {
                     ret4 <-  read_csv(inFile4$datapath, col_types = readr::cols())
                     d <- 1
                     setProgress(1, detail = paste("Completed"))
                   }
      )
    }
    
    if (!is.null(inFile5))  {
      withProgress(message = 'Reading the data',
                   detail = "This may take a while", value = 0, {
                     ret5 <-  read_csv(inFile5$datapath, col_types = readr::cols())
                     e <- 1
                     setProgress(1, detail = paste("Completed"))
                   }
      )
    }
    
    
    if (a + b + c + d + e < 2) {
      sendSweetAlert(
        session = session,
        title =  "Missing DEGdata",
        text = paste0("Please upload at least 2 DEG files"),
        type = "error"
      )
      return(NULL)
    }
    return(list("DEG1" = ret1,
                "DEG2" = ret2,
                "DEG3" = ret3,
                "DEG4" = ret4,
                "DEG5" = ret5))
  })  
  getFormula <- reactive({
    form <- NULL
    cond <- input$condition
    cov <- input$covariates
    if(is.null(cov)){
      if(str_length(cond) > 0) form <- as.formula(paste0("~ ", cond))
    } else if(str_length(cond) > 0 & (length(cov) == 0)) {
      form <- as.formula(paste0("~ ", cond))
    } else if(length(cov) > 0) {
      form <- as.formula(paste0("~ ",paste(cov,collapse = "+")," + ", cond))
    }
    return(form)
  })
  
  output$formulatext <- renderText({
    f <- getFormula()
    if(!is.null(f)) return(as.character(f))
    return("")
  })
  
  output$metadata.tbl <-  DT::renderDataTable({
    metadata <- readMetaData()
    if(is.null(metadata)) {
      return(NULL)
    }
    metadata  %>% createTable2(show.rownames=F)
  })
  
  observeEvent(input$condition, {
    metadata <- readMetaData()
    if(!is.null(metadata)) {
      updateSelectizeInput(session, 'reference', choices =  as.character(unique(metadata %>% pull(input$condition))), server = TRUE)
    }
  })
  
  observeEvent(input$dea, {
    updateTabsetPanel(session, inputId="DEA", selected = "DEA results")
    if(!is.null(get.DEA.results()) & input$DEmethod=="deseq") updateSelectizeInput(session, 'deaSelect', 
                                                         choices =  resultsNames(get.DEA.results())[-1],
                                                         selected = resultsNames(get.DEA.results())[-1][1], 
                                                         server = TRUE)
    if(!is.null(get.DEA.results()) & input$DEmethod=="edger") updateSelectizeInput(session, 'deaSelect', 
                                                                                   choices =  dimnames(get.DEA.results())[[2]],
                                                                                   selected = dimnames(get.DEA.results())[[2]][1], 
                                                                                   server = TRUE)
    output$dea.results <-  DT::renderDataTable({
      res <- get.DEA.results()
      if(is.null(res)) return(NULL)
      deaSelect <- input$deaSelect
      if (input$DEmethod =="deseq") { 
        if(str_length(deaSelect) == 0) {
          if(input$log2FoldChange > 0){
            tbl <-  as.data.frame(results(res,
                                          lfcThreshold = input$log2FoldChange,  
                                          altHypothesis = "greaterAbs"))
            
          } else{
            tbl <-  as.data.frame(results(res,name = deaSelect))
          }
        } else {
          if(input$lfc) {
            tbl <-  as.data.frame(lfcShrink(res, coef = deaSelect)) ### adding apeglm option breaks results, remove to make sure
            
            withProgress(message = 'Shrink log2 fold changes',
                         detail = "Shrinking...", value = 0, {
                           tbl <-  as.data.frame(lfcShrink(res, 
                                                           type = "apeglm", 
                                                           coef = deaSelect,
                                                           svalue=TRUE)) ### adding apeglm option breaks results
                         })
            
          } else {
            if(input$log2FoldChange > 0){
              tbl <-  as.data.frame(results(res,
                                            name = deaSelect,
                                            lfcThreshold = input$log2FoldChange,  
                                            altHypothesis = "greaterAbs"))
            } else {
              tbl <-  as.data.frame(results(res,name = deaSelect))
            }
          }
        }
        tbl <- cbind("Symbol" = rownames(tbl), tbl)
        tbl %>% createTable2(show.rownames = F)}
      #Data processing for edgeR method
      else if (input$DEmethod =="edger") {
        if(str_length(deaSelect) == 0) {
          if(input$log2FoldChange > 0){
            tbl <-  as.data.frame(res$table)
            
          } else{
            tbl <-  as.data.frame(res$table)
          }
        } else {
          if(input$lfc) {
            tbl <-  as.data.frame(res$table)
          } else {
            if(input$log2FoldChange > 0){
              tbl <-  as.data.frame(res$table)
            } else {
              tbl <-  as.data.frame(res$table)
            }
          }
        }
        tbl <- cbind("Symbol" = rownames(tbl), tbl)
        tbl %>% createTable2(show.rownames = F)}        
      
    })
  })
  
  #
  # Pathway analysis tab
  #  source: https://guangchuangyu.github.io/pathway-analysis-workshop/
  #  
  
  # Return list of DEA genes sorted and the names of the most significant ones  
  readDEA <- reactive({
    data <- NULL
    inFile <- input$deafile
    gmtFile <- input$gmtFile
    if (!is.null(inFile))  {
      withProgress(message = 'Reading the data',
                   detail = "This may take a while", value = 0, {
                     data <-  read_csv(inFile$datapath, col_types = readr::cols())
                     setProgress(1, detail = paste("Completed"))
                   }
      )
      if(!is.data.frame(data)){
        withProgress(message = 'Reading the data',
                     detail = "This may take a while", value = 0, {
                       data <-  read_csv2(inFile$datapath, col_types = readr::cols())
                       setProgress(1, detail = paste("Completed"))
                     }
        )
      }
    }
    if (!is.null(gmtFile)) {
      withProgress(message = 'Reading the GMT data',
                   detail = "This may take a while", value = 0, {
                     pathway <- gmtPathways(gmtFile$datapath)
                     gmtDF <- read_tsv(gmtFile$datapath,
                                       col_names = c("pathway", "description"),
                                       col_types = readr::cols())
                     gmtDF <- gmtDF[,1:2]
                     setProgress(1, detail = paste("Completed"))
                   }
      )
    }
    if(is.null(data)) {
      sendSweetAlert(
        session = session,
        title = "Missing input data", 
        text = "Please upload DEA results",
        type = "error"
      )
      return(NULL)
    }
    
    
    if(!"log2FoldChange" %in% colnames(data)){
      sendSweetAlert(
        session = session,
        title = "Data input not as expected", 
        text = paste0("No log2FoldChange column in the input"),
        type = "error"
      )
      return(NULL)
    } 
    
    
    if(!"pvalue" %in% colnames(data) & !"svalue" %in% colnames(data)){
      sendSweetAlert(
        session = session,
        title = "Data input not as expected", 
        text = paste0("No pvalue or svalue column in the input"),
        type = "error"
      )
      return(NULL)
    } 
    
    ret <- subset(data, log2FoldChange>=2.0 | log2FoldChange <= -2.0 )
   
    if("svalue" %in% colnames(ret)) {
      ret <- subset(ret, svalue <= 0.05)
    }
    else if("pvalue" %in% colnames(ret)) {
      ret <- subset(ret, pvalue <= 0.05)
    }

    # For GSEA
    
    if(input$earankingmethod == "log Fold Change") {
      geneList.metric <- ret$log2FoldChange
    } else if(input$earankingmethod ==  "-log10(P-value) * sig(log2FC)") {
      if ("pvalue" %in% colnames(data)) {
        geneList.metric <- -log10(ret$pvalue) * sign(ret$log2FoldChange)        
      }
      else {
        geneList.metric <- -log10(ret$svalue) * sign(ret$log2FoldChange)         
      }

    } else {
      if ("pvalue" %in% colnames(data)) {
        geneList.metric <- -log10(ret$pvalue) * ret$log2FoldChange        
      }
      else {
        geneList.metric <- -log10(ret$svalue) * ret$log2FoldChange        
      }

    }
    geneList.metric.ensembl <- geneList.metric
    names(geneList.metric.ensembl) <- ret$Symbol
    geneList.metric.ensembl <- sort(geneList.metric.ensembl, decreasing = TRUE)
    geneList.metric.ensembl <- geneList.metric.ensembl[!is.na(names(geneList.metric.ensembl))]
    
    names(geneList.metric) <- ret$Symbol
    geneList.metric <- sort(geneList.metric, decreasing = TRUE)
    geneList.metric <- geneList.metric[!is.na(names(geneList.metric))]
    geneList.metric <- geneList.metric[!duplicated(names(geneList.metric))]
    return(list("geneList" = geneList.metric,
                "geneList.ensembl" = geneList.metric.ensembl,
                "pathway" = pathway,
                "dea.results" = ret,
                "gmtDF" = gmtDF))
  })

  # Perform selected analysis
  # 
  observeEvent(input$deaanalysisselect,{
    if(input$deaanalysisselect ==  "MSigDb analysis") {
      shinyjs::show(id = "msigdbtype", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    } else {
      shinyjs::hide(id = "msigdbtype", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    }    
    
    if(input$deaanalysisselect ==  "Gene Ontology Analysis") {
      shinyjs::show(id = "gotype", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    } else {
      shinyjs::hide(id = "gotype", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    }  
  })
  
  observeEvent(input$deaanalysistype,  {
    if(input$deaanalysistype == "GSEA") {
      updateSelectizeInput(session, 'ea_plottype', 
                           selected = "Dot plot",
                           choices = c("Dot plot",
                                       "GSEA table plot"),
                           server = TRUE)
      shinyjs::show(id = "eagsearankingui", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
      shinyjs::hide(id = "eaorasectui", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    } else {
      updateSelectizeInput(session, 'ea_plottype', 
                           selected = "Dot plot",
                           choices = c("gost plot",
                                       "GSEA table plot",
                                       "Running score and preranked list",
                                       "Ranked list of genes",
                                       "Enrichment map (network)"),
                           server = TRUE)
      shinyjs::show(id = "eaorasectui", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
      shinyjs::hide(id = "eagsearankingui", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
      
    }
  })
  
  observeEvent(input$enrichementbt,  {
    enrichement.analysis <- reactive({
      data <- readDEA()
      if(is.null(data)) {
        
        return(NULL)
      }
      
#      if(length(data$dea.genes) == 0){
        
#        sendSweetAlert(
#          session = session,
#          title = "No genes identified", 
#          text = "We could not map the genes Symbols to entrez gene ID. Please check input data.",
#          type = "error"
#        )
        
#        return(NULL)
#      }
      if(isolate({input$deaanalysisselect}) != "Gene Ontology Analysis") {
#        dea.genes <- data$dea.genes
        geneList <- data$geneList
      } else {
#        dea.genes <- data$dea.genes.ensembl
        geneList <- data$geneList.ensembl
        
      }
      withProgress(message = 'Performing analysis',
                   detail = "It might take a while...", value = 0, {
                     
                     if(isolate({input$deaanalysisselect}) == "Gene Ontology Analysis"){
                       message("o Gene Ontology Analysis")
                       
                       if (isolate({input$deaanalysistype})=="gprofiler") {
                         print(type(names(geneList)))
                         results <- gost(names(geneList), organism="ncrassa")
                         print(results)
                       }
                       
                       else {
                         results <- fgsea(data$pathway, data$geneList, 
                                          minSize=4,
                                          maxSize=500, 
                                          nperm=input$permutations,
                                          gseaParam=input$gseaParam) %>%
                           as_tibble() %>%
                           arrange(padj)
                         results <- results %>%
                           as_tibble() %>%
                           arrange(desc(NES))
                         results %>%
                           dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
                           arrange(pval) %>%
                           DT::datatable()
                         results <- merge(results, data$gmtDF, by="pathway")
                       }
                       
                     } else if(isolate({input$deaanalysisselect}) == "KEGG Analysis"){
                       message("o KEGG Analysis")
                       if(isolate({input$deaanalysistype}) == "ORA"){
                         results <- enrichKEGG(dea.genes, 
                                               pvalueCutoff = input$enrichmentfdr,
                                               organism = "hsa")
                       } else {
                         results <- gseKEGG(geneList, 
                                            pvalueCutoff = input$enrichmentfdr,
                                            organism = "hsa", 
                                            nPerm = 10000)
                       }
                     } 
                   })
      return(results)
    })  
    output$tbl.analysis <-  DT::renderDataTable({
      input$enrichementbt
      tbl <- enrichement.analysis()
      if(is.null(tbl)) {
        
        return(NULL)
      }
      tbl %>% createTable2(show.rownames = F)
#summary %>%       
    })

    
    updateSelectizeInput(session, 'gsea_gene_sets', 
                         selected = enrichement.analysis()$Description[1],
                         choices =  enrichement.analysis()$Description,
                         server = TRUE)
    
    getEnrichementPlot <- reactive({
      results <- enrichement.analysis()
      if(is.null(results)) return(NULL)
      p <- NULL
      
      if(nrow(summary(results)) == 0){
        if("pvalueCutoff"  %in% slotNames(results)){
          aux <- results@pvalueCutoff
        } else {
          aux <- results@params$pvalueCutoff
        }
        
        sendSweetAlert(
          session = session,
          title =  "No enriched terms found", 
          text =  paste0("No results for enrichment analysis P-value cut-off = ", aux),
          type = "error"
        )
        return(NULL)
      }
      
      
      if(isolate({input$deaanalysistype}) == "GSEA") {
        # GSEA plots
        
        if( input$ea_plottype == "Dot plot") {
          resCutoff <- isolate({input$enrichmentfdr})
          res <- subset(results, pval<resCutoff)
          legLabel <- paste("adj.p < ", resCutoff,sep="")
          p <- ggplot(res, aes(reorder(description, NES), NES)) +
            geom_col(aes(fill=padj<isolate({input$enrichmentfdr}))) +
            coord_flip() +
            labs(x="Pathway", y="Normalized Enrichment Score",
                 title="Hallmark pathways NES from GSEA",
                 fill=legLabel) + 
            theme_minimal()
        }
        else if(input$ea_plottype == "GSEA table plot") {
          data <- readDEA()
          resCutoff <- isolate({input$enrichmentfdr})
          res <- results
          topPathways <- res %>% 
            top_n(input$ea_nb_categories, wt=-padj) %>% 
            arrange(-NES) %>% 
            pull(pathway)
          p <- plotGseaTable(data$pathway[topPathways], 
                        data$geneList, 
                        res, 
                        gseaParam = input$gseaParam)
        }

        if(input$ea_plottype == "GSEA Plot") {
          p <- ridgeplot(results,showCategory = input$ea_nb_categories)
        }
        
        if(input$ea_plottype ==  "Running score and preranked list") {
          p <- gseaplot2(results, geneSetID = match(input$gsea_gene_sets, results$Description))
        }
        if(input$ea_plottype ==  "Ranked list of genes") {
          p <- gsearank(results,  which(input$gsea_gene_sets ==  results$Description), 
                        title = results[ match(input$gsea_gene_sets, results$Description), "Description"])
          
          pp <- lapply( match(input$gsea_gene_sets, results$Description), function(i) {
            anno <- results[i, c("NES",  "p.adjust")]
            lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
            
            es <- results[i, "enrichmentScore"]
            x <- ifelse(es < 0, 0, length(geneList) * .8)
            gsearank(results, i, results[i, 2]) + xlab(NULL) +ylab(NULL) +
              annotate("text", x,  es * .5, label = lab, 
                       hjust=0, vjust=0, size = 4) + xlim(0, 12500)
          })
          p <- plot_grid(plotlist=pp, ncol=1)
        }
        if( input$ea_plottype == "Enrichment map (network)") {
          p <- emapplot(results)
          
        }          
      }
      else if(isolate({input$deaanalysistype})=="gprofiler") {
        if(input$ea_plottype == "Dot plot") {
          p <- gostplot(results) 
        }
      }
      p
    })
    
    # Save figure 
    output$saveenrichementpicture <- downloadHandler(
      filename = function(){input$enrichementPlot.filename},
      content = function(file) {
        if(tools::file_ext(input$enrichementPlot.filename) == "png") {
          device <- function(..., width, height) {
            grDevices::png(..., 
                           width = isolate({input$ea_width}), 
                           height = isolate({input$ea_height}), 
                           res = 300, 
                           units = "in")
          }
        } else if(tools::file_ext(input$enrichementPlot.filename) == "pdf") {
          device <- "pdf"
        } else if(tools::file_ext(input$enrichementPlot.filename) == "svg") {
          device <- function(..., width, height) {
            grDevices::svg(..., width = isolate({input$ea_width}), height = isolate({input$ea_height}))
          } 
        } else {
          
          sendSweetAlert(
            session = session,
            title = "Extension not recognized", 
            text =  paste0("svg, pdf and png allowed"),
            type = "error"
          )
        }
        p <- getEnrichementPlot()
        ggsave(file, 
               plot = p , 
               device = device, 
               width = isolate({input$ea_width}), 
               height = isolate({input$ea_height}), 
               units = "in")
      })
    
    output$plotenrichment <- renderPlot({
      getEnrichementPlot()
    })
    
    
    
  })
  
  
  output$downloadExampleDEAData <- downloadHandler(
    filename = function() {
      "subtype_BRCA_Subtype_PAM50_LumA_vs_Basal.csv"
    },
    content = function(file) {
      metadata <- readr::read_csv("test/subtype_BRCA_Subtype_PAM50_LumA_vs_Basal.csv")
      write_csv(metadata, file)
    }
  )
  
  
  
  
  observeEvent(input$ea_plottype, {
    if(input$ea_plottype %in% c("Running score and preranked list", "Ranked list of genes")){
      shinyjs::show(id = "gsea_gene_sets", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
      shinyjs::hide(id = "ea_nb_categories", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    } else {
      shinyjs::hide(id = "gsea_gene_sets", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
      shinyjs::show(id = "ea_nb_categories", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    }
  })
  
  #---------------------------------
  # Volcano plot tab
  #---------------------------------
  output$volcanoplot <- renderPlotly({
    res <- get.DEA.results()
    if(is.null(res)) return(NULL)
    
    if(input$DEmethod=="deseq") {
      x.cut <- input$log2FoldChange
      y.cut <- input$padj
      
      deaSelect <- input$deaSelect
      if(str_length(deaSelect) == 0) {
        dea <-  as.data.frame(results(res))
      } else {
        if(input$lfc) {
          
          #dea <-  as.data.frame(lfcShrink(res, coef = deaSelect)) ### comment out this line to merge with tiagos changes, adding in apeglm for lfcshrink, doesn't break the app
          
          withProgress(message = 'Shrink log2 fold changes',
                       detail = "Shrinking...", value = 0, {
                         dea <-  as.data.frame(lfcShrink(res, 
                                                         coef = deaSelect,
                                                         type = "apeglm",
                                                         svalue=TRUE)) ### adding apeglm option breaks results
                       }
          )
          
        } else {
          if(x.cut > 0){
            dea <-  as.data.frame(results(res,
                                          name = deaSelect,
                                          lfcThreshold = x.cut,
                                          altHypothesis = "greaterAbs"))
          } else {
            dea <-  as.data.frame(results(res, name = deaSelect))
          }
        }
      }
      
      dea$group <- "Not Significant"
      dea[which(dea$padj < y.cut & dea$log2FoldChange < -x.cut ),"group"] <- "Downregulated"
      dea[which(dea$padj < y.cut & dea$log2FoldChange > x.cut ),"group"] <- "Upregulated"
      f <- list(
        family = "Courier New, monospace",
        size = 18,
        color = "#7f7f7f"
      )
      x <- list(
        title = "log2FoldChange",
        titlefont = f
      )
      y <- list(
        title = "-log10(p-value adjusted)",
        titlefont = f
      )
      
      if (input$lfc) {
        dea[which(dea$svalue < y.cut & dea$log2FoldChange < -x.cut ),"group"] <- "Downregulated"
        dea[which(dea$svalue < y.cut & dea$log2FoldChange > x.cut ),"group"] <- "Upregulated"
        p <- plot_ly(data = dea, 
                     x = dea$log2FoldChange, 
                     y = -log10(dea$svalue), 
                     text = rownames(dea), 
                     mode = "markers", 
                     color = dea$group) %>% 
          layout(title ="Volcano Plot") %>%
          layout(xaxis = x, yaxis = y)  %>%
          layout(shapes=list(list(type='line', 
                                  x0 = x.cut, 
                                  x1 = x.cut, 
                                  y0 = 0, 
                                  y1 = max(-log10(dea$svalue),na.rm = T), 
                                  line=list(dash='dot', width=1)),
                             list(type='line', 
                                  x0 = -x.cut, 
                                  x1 = -x.cut, 
                                  y0 = 0, 
                                  y1 = max(-log10(dea$svalue),na.rm = T), 
                                  line =list(dash='dot', width=1)),
                             list(type ='line', 
                                  x0 = min(dea$log2FoldChange), 
                                  x1 = max(dea$log2FoldChange), 
                                  y0 =  -log10(y.cut), 
                                  y1 =  -log10(y.cut), 
                                  line = list(dash='dot', width=1))
          ) 
          )        
      }
      else {
        dea[which(dea$padj < y.cut & dea$log2FoldChange < -x.cut ),"group"] <- "Downregulated"
        dea[which(dea$padj < y.cut & dea$log2FoldChange > x.cut ),"group"] <- "Upregulated"
        p <- plot_ly(data = dea, 
                     x = dea$log2FoldChange, 
                     y = -log10(dea$padj), 
                     text = rownames(dea), 
                     mode = "markers", 
                     color = dea$group) %>% 
          layout(title ="Volcano Plot") %>%
          layout(xaxis = x, yaxis = y)  %>%
          layout(shapes=list(list(type='line', 
                                  x0 = x.cut, 
                                  x1 = x.cut, 
                                  y0 = 0, 
                                  y1 = max(-log10(dea$padj),na.rm = T), 
                                  line=list(dash='dot', width=1)),
                             list(type='line', 
                                  x0 = -x.cut, 
                                  x1 = -x.cut, 
                                  y0 = 0, 
                                  y1 = max(-log10(dea$padj),na.rm = T), 
                                  line =list(dash='dot', width=1)),
                             list(type ='line', 
                                  x0 = min(dea$log2FoldChange), 
                                  x1 = max(dea$log2FoldChange), 
                                  y0 =  -log10(y.cut), 
                                  y1 =  -log10(y.cut), 
                                  line = list(dash='dot', width=1))
          ) 
          )
    }
  }
    if(input$DEmethod=="edger") {
      x.cut <- input$log2FoldChange
      y.cut <- input$padj
      
      deaSelect <- input$deaSelect
      if(str_length(deaSelect) == 0) {
        dea <-  res$table
      } else {
          if(x.cut > 0){
            dea <-  select(res$table, c(matches(deaSelect),"F","PValue"))
          } else {
            dea <-  select(res$table, c(matches(deaSelect),"F","PValue"))
          }
        
      }
      names(dea)[names(dea) == deaSelect] <- "logFC"
      dea$group <- "Not Significant"
      dea[which(dea$PValue < y.cut & dea$logFC < -x.cut ),"group"] <- "Downregulated"
      dea[which(dea$PValue < y.cut & dea$logFC > x.cut ),"group"] <- "Upregulated"
      
      
      f <- list(
        family = "Courier New, monospace",
        size = 18,
        color = "#7f7f7f"
      )
      x <- list(
        title = "log2FoldChange",
        titlefont = f
      )
      y <- list(
        title = "-log10(p-value adjusted)",
        titlefont = f
      )
      p <- plot_ly(data = dea, 
                   x = dea$logFC,
                   y = -log10(dea$PValue), 
                   text = rownames(dea), 
                   mode = "markers", 
                   color = dea$group) %>% 
        layout(title ="Volcano Plot") %>%
        layout(xaxis = x, yaxis = y)  %>%
        layout(shapes=list(list(type='line', 
                                x0 = x.cut, 
                                x1 = x.cut, 
                                y0 = 0, 
                                y1 = max(-log10(dea$PValue),na.rm = T), 
                                line=list(dash='dot', width=1)),
                           list(type='line', 
                                x0 = -x.cut, 
                                x1 = -x.cut, 
                                y0 = 0, 
                                y1 = max(-log10(dea$PValue),na.rm = T), 
                                line =list(dash='dot', width=1)),
                           list(type ='line', 
                                x0 = min(dea$logFC), 
                                x1 = max(dea$logFC), 
                                y0 =  -log10(y.cut), 
                                y1 =  -log10(y.cut), 
                                line = list(dash='dot', width=1))
        ) 
        )
    }
    return(p)
  })
  #---------------------------------
  # Venn Diagram plot tab
  #---------------------------------
  output$vennplot <- renderPlotly({
    input$venn
    log2Cutoff <- input$vennLog2FoldChange
    pvalCutoff <- input$vennPadj
    data <- get.DEG.results()
    DEG1 <- NULL
    DEG2 <- NULL
    DEG3 <- NULL
    DEG4 <- NULL
    if(is.null(data)) return(NULL)
    if(!is.null(data$DEG1)) {
      DEG1 <- subset(data$DEG1, svalue<pvalCutoff,
                          select=c(Symbol,log2FoldChange))
      DEG1 <- subset(DEG1, log2FoldChange > log2Cutoff | log2FoldChange < -log2Cutoff,
                     select=Symbol)
      DEG1$Group <- input$DEfile1Name
    }
    if(!is.null(data$DEG2)) {
      DEG2 <- subset(data$DEG2, svalue<pvalCutoff,
                     select=c(Symbol,log2FoldChange))
      DEG2 <- subset(DEG2, log2FoldChange > log2Cutoff | log2FoldChange < -log2Cutoff,
                     select=Symbol)
      DEG2$Group <- input$DEfile2Name
    }
    if(!is.null(data$DEG3)) {
      DEG3 <- subset(data$DEG3, svalue<pvalCutoff,
                     select=c(Symbol,log2FoldChange))
      DEG3 <- subset(DEG3, log2FoldChange > log2Cutoff | log2FoldChange < -log2Cutoff,
                     select=Symbol)
      DEG3$Group <- input$DEfile3Name
    }
    if(!is.null(data$DEG4)) {
      DEG4 <- subset(data$DEG4, svalue<pvalCutoff,
                     select=c(Symbol,log2FoldChange))
      DEG4 <- subset(DEG4, log2FoldChange > log2Cutoff | log2FoldChange < -log2Cutoff,
                     select=Symbol)
      DEG4$Group <- input$DEfile4Name
    }
    allData <- rbind(DEG1, DEG2, DEG3,DEG4)
    allData <- drop_na(allData)
    vennCat <- unique(allData$Group)

###################################################################
    p <- plotVenn(vennCat, category=vennCat, lty="blank", fill = c("skyblue", "pink1"))
    return(p)
  })

  
  
}