source("library.R")$value

ui <- fluidPage(title = "RNA-UI",
                theme = shinytheme("superhero"),
                tags$head(tags$style(
                  HTML('
                       #sidebar {
                       background-color: #FFFFFF;
                       }

                       body, label, input, button, select { 
                       font-family: "Arial";
                       }
                       
                       .sw-dropdown-in {
                       background-color: #2b3e50;
                       margin: 0;
                       }

                       .btn-file {
                        background-color:#5B81AE; 
                        border-color: #5B81AE; 
                        background: #5B81AE;
                       }
                      
         
                      .bttn-bordered.bttn-sm {
                          width: 200px;
                          text-align: left;
                          margin-bottom : 0px;
                          margin-top : 20px;
                       }
                       '
                  )
                )),
                titlePanel("RNA-UI"),
                useShinyjs(),
                tabsetPanel( #type = "pills",
                  tabPanel("Upload Data",  
                           icon = icon("table"),
                           column(2,
                                  #sidebarPanel(id="sidebar",width = 3,
                                  dropdown(label = "Data upload",
                                           icon = icon("file-csv"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           fileInput("rawcounts", "Choose CSV File",
                                                     multiple = TRUE,
                                                     accept = c("text/csv",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv")),
                                           tags$div(
                                             HTML(paste(help_text))
                                           )
                                  ),
                                  
                                  
                                  #Changed from 'Table Selection'
                                  dropdown(label = "Data Normalization",
                                           icon = icon("table"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           selectInput("select_tab1", "Select Transform", transforms, multiple = FALSE)
                                  ), ##need individual selectInputs for each tab
                                  dropdown(label = "Download data",
                                           icon = icon("download"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           downloadButton("downloadNormalizedData", "Download normalized files",class = "btn-primary")
                                  ),
                                  dropdown(label = "Generate report",
                                           icon = icon("file-code"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           downloadButton("reportNorm", "Download report",class = "btn-primary")
                                  ),
                                  
                                  
                                  dropdown(label = "Gene selection",
                                           icon = icon("mouse-pointer"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           fileInput("input_gene_list_tab1", "Input Gene Symbol List (Optional)", 
                                                     multiple = FALSE, 
                                                     accept = NULL, 
                                                     width = NULL, 
                                                     buttonLabel = "Browse", 
                                                     placeholder = "No file selected"), ##how to increase max upload size
                                           textAreaInput(inputId = "input_gene_list_area",
                                                         label = "Gene list filter: separate gene names by , or ; or newline",
                                                         value =  "", 
                                                         width = "100%"),
                                           actionButton("input_gene_list_but", 
                                                        "Select Rows",
                                                        width = "100%", 
                                                        class = "btn-primary"), ##do this to put selected rows at top of data table, trying it out
                                           actionButton("select_most_variable", 
                                                        "Select 1000 genes of highest variance",
                                                        width = "100%", 
                                                        class = "btn-primary"), ##do this to put selected rows at top of data table, trying it out
                                           actionButton("unselect_all", 
                                                        "Deselect all genes",
                                                        width = "100%", 
                                                        class = "btn-primary")
                                  ), 

                                  dropdown(label = "Sample selection",
                                           icon = icon("mouse-pointer"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           fileInput("input_sample_list_tab1", "Input Sample List (Optional)", 
                                                     multiple = FALSE, 
                                                     accept = NULL, 
                                                     width = NULL, 
                                                     buttonLabel = "Browse", 
                                                     placeholder = "No file selected"),
                                           
                                           selectInput("selectedSamples", "Select Samples for Comparison", NULL, 
                                                       multiple = TRUE)
                                  )
                           ),
                           column(10,
                                  DT::dataTableOutput('tbl.tab1') 
                           )
                  ),
                  
                  #Probably need to remove this. We will see ###################################################
                  tabPanel("Visualization", ##changing from tab 2, but still usibg tab2 in other parts of code
                           #icon = icon("object-group"),
                           icon = icon("image"),
                           tabsetPanel(type = "pills",
                                       tabPanel("Expression plots",
                                                icon = icon("bar-chart-o"),
                                                bsAlert("genemessage"),
                                                column(2,
                                                       selectInput("Expheatcolor", "Color Sample Treatments by:", NULL, multiple = FALSE)),
                                                hidden(
                                                  div(id = "expression_plots",
                                                      h3('Expression Barplot'), 
                                                      plotlyOutput("barplot", width = "100%") %>% withSpinner(type = 6)
                                                  )),
                                                hidden(
                                                  div(id = "expression_heatmap",
                                                      h3('Expression Heatmap'),
                                                      #selectInput("select_z_score", 
                                                      #            label = "Standardized scores?", 
                                                      #            choices = c("No","Rows z-score", "Columns z-score"), 
                                                      #            multiple = FALSE),
                                                      iheatmaprOutput("heatmap_expr",height = "auto") %>% withSpinner(type = 6)
                                                  )
                                                )
                                       ),
                                       tabPanel("Clustering plots",
                                                icon = icon("object-group"),
                                                div(id = "cluster_plots",
                                                    column(2,
                                                           h3('Correlation Heatmap'), 
                                                           selectInput("select_clus_type", 
                                                                       label = "Cluster correlation", 
                                                                       choices = c("Sample","Genes"), 
                                                                       multiple = FALSE),
                                                           selectInput("select_clus", "Cluster by what genes", 
                                                                       c("All genes", "Selected genes"), 
                                                                       multiple = FALSE)
                                                    ),
                                                    column(9,
                                                           iheatmaprOutput("heatmap_clus",height = "800px") %>% withSpinner(type = 6)
                                                    )
                                                )
                                       ), tabPanel("PCA plots",
                                                   icon = icon("object-group"),
                                                   div(id = "pca_plots",
                                                       bsAlert("genemessage3"),
                                                       column(2,
                                                              selectInput("select_pca_type", 
                                                                          label = "PCA genes", 
                                                                          choices = c("Top 1000 variable genes", "All genes", "Selected genes"), 
                                                                          multiple = FALSE),
                                                              selectInput("pca_dimensions", 
                                                                          label = "Number of dimensions", 
                                                                          choices = c("2D", "3D"), 
                                                                          multiple = FALSE),
                                                              selectInput("pcacolor", "Color samples by", NULL, multiple = FALSE),
                                                              radioButtons("confEllipse", NULL,
                                                                          choices=c("No Confidence Interval"= "noConf",
                                                                                    "95% Confidence Interval" = "confEll")),
                                                              radioButtons("geneLoadings", NULL,
                                                                           choices=c("No Gene Vectors"= "noGeneLoad",
                                                                                     "Include Gene Loading Vectors" = "GeneLoad")),
                                                              sliderInput("topGenes", "Top Contributing Genes to Include",
                                                                          min=1, max=1000, value=20),
                                                              downloadButton("downloadPCA", "Download Principal Components", class = "btn-primary")),
                                                       column(6,
                                                              plotlyOutput("pca_plot",height = "600",width = "600") %>% withSpinner(type = 6)
                                                       )

                                                   ),
                                       )
                           )
                  ),
                  tabPanel("Differential Expression Analysis", 
                           icon = icon("calculator"),
                           column(2,
                                  #sidebarPanel(id="sidebar",width = 3,
                                  dropdown(label = "Metadata upload",
                                           icon = icon("file-csv"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeft,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           # Input: Select a file ----
                                           downloadButton('downloadData', 'Download example metadata file',class = "btn-primary"),
                                           fileInput("metadata", "Choose CSV File",
                                                     multiple = TRUE,
                                                     accept = c("text/csv",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv")),
                                           tags$div(
                                             HTML(paste(help_text2))
                                           )
                                  ),
                                  dropdown(label = "DE analysis",
                                           icon = icon("calculator"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           
                                           selectInput("condition", "Select condition column for DEA", NULL, multiple = FALSE), ##need individual selectInputs for each tab
                                           fileInput("input_sample_list_DEA", "Input Sample List (Optional)", 
                                                     multiple = FALSE, 
                                                     accept = NULL, 
                                                     width = NULL, 
                                                     buttonLabel = "Browse", 
                                                     placeholder = "No file selected"),
                                           selectInput("covariates", 
                                                       label = "Select covariates for DEA",
                                                       choices =  NULL, 
                                                       multiple = TRUE), ##need individual selectInputs for each tab
                                           verbatimTextOutput("formulatext"),
                                           selectInput("reference", "Select reference level for DEA", NULL, multiple = FALSE), ##need individual selectInputs for each tab
                                           selectInput("DEmethod", "Select DE Method", 
                                                       choices= c("DESeq2"="deseq",
                                                                  "edgeR"="edger"), 
                                                       multiple=FALSE),
                                           checkboxInput(inputId="lfc", label = "apeGLM NOTE: only perform on DESeq", value = FALSE, width = NULL),
                                           actionButton("dea", "Perform DEA", class = "btn-primary")
                                  ),
                                  dropdown(label = "Select Results",
                                           icon = icon("table"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           selectInput("deaSelect", "Select results", NULL, multiple = FALSE), ##need individual selectInputs for each tab
                                           downloadButton("downloadDEAFiles", "Download DEA Results",class = "btn-primary")
                                  ),
                                  dropdown(label = "Volcano plot",
                                           icon = icon("mountain"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           numericInput("log2FoldChange", "log2FoldChange  cut-off:", value = 0, min = 0, max = 10, step = 0.1),
                                           numericInput("padj", "P adjusted cut-off or s-value:", 0.01, min = 0, max = 1,step = 0.1)
                                  ),
                                  dropdown(label = "Generate report",
                                           icon = icon("file-code"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           downloadButton("reportDEA", "Download Report",class = "btn-primary")
                                  )
                           ),
                           column(10,
                                  tabsetPanel(type = "pills",
                                              id = "DEA",
                                              tabPanel("Metadata",
                                                       tags$hr(),
                                                       DT::dataTableOutput('metadata.tbl')
                                              ), 
                                              tabPanel("DEA results",
                                                       tags$hr(),
                                                       DT::dataTableOutput('dea.results') 
                                              ),
                                              tabPanel("Volcano plot",
                                                       tags$hr(),
                                                       plotlyOutput('volcanoplot') %>% withSpinner(type = 6)
                                              )
                                  )
                           )
                  ),
                  tabPanel("Enrichment analysis", 
                           icon = icon("calculator"),
                           column(2,
                                  #sidebarPanel(id="sidebar",width = 3,
                                  dropdown(label = "DESeq results upload      ",
                                           icon = icon("file-csv"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           downloadButton('downloadExampleDEAData', 'Download example DEA file',class = "btn-primary"),
                                           fileInput("deafile", "Choose CSV File",
                                                     multiple = TRUE,
                                                     accept = c("text/csv",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv"))
                                  ),
                                  dropdown(label = "GMT file upload      ",
                                           icon = icon("file-import"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           downloadButton('downloadExampleGMT', 'Download example GMT file',class = "btn-primary"),
                                           fileInput("gmtFile", "Choose GMT File",
                                                     multiple = TRUE,
                                                     accept = c("text/tab",
                                                                "text/tab-separated-values,text/plain",
                                                                ".gmt"))
                                  ),
                                  
                                  #Need to add different options for the CLUEGO, pathway enrichment, settings,
                                  #also need to be able to upload custom GMT files and select necessary enrichment
                                  #files for input into CLUEGO
                                  dropdown(label = "Enrichment Analysis",
                                           icon = icon("calculator"),
                                           style = "bordered", 
                                           size =  "sm",
                                           status = "primary", 
                                           width = "300px",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           selectInput("deaanalysistype", 
                                                       "Select the type of analysis", 
                                                       c("GSEA (gene set enrichment analysis)" = "GSEA",
                                                         "g:Profiler2 (over enrichment analysis)" = "gprofiler"), 
                                                       multiple = FALSE),
                                           selectInput("deaanalysisselect", 
                                                       "Select the analysis", 
                                                       c("Gene Ontology Analysis",
                                                         "KEGG Analysis"), 
                                                       multiple = FALSE),
                                           selectInput("msigdbtype", 
                                                       "Select collection for Molecular Signatures Database", 
                                                       c("All human gene sets" = "All",
                                                         "H: hallmark gene sets" = "H",
                                                         "C1: positional gene sets" = "C1",
                                                         "C2: curated gene sets" = "C2",
                                                         "C3: motif gene sets" = "C3",
                                                         "C4: computational gene sets" = "C4",
                                                         "C5: GO gene sets" = "C5",
                                                         "C6: oncogenic signatures" = "C6",
                                                         "C7: immunologic signatures" = "C7"), 
                                                       multiple = FALSE),
                                           selectInput("gotype", 
                                                       "Select collection for Molecular Signatures Database", 
                                                       c("Molecular Function"="MF",
                                                         "Cellular Component"="CC",
                                                         "Biological Process" = "BP"), 
                                                       multiple = FALSE),
                                           numericInput("enrichmentfdr", 
                                                        "P-value cut-off:", 
                                                        value = 0.05, 
                                                        min = 0, 
                                                        max = 1, 
                                                        step = 0.05),
                                           numericInput("permutations", 
                                                        "Number of permutations (rec. 1000)", 
                                                        value = 1000, 
                                                        min = 100, 
                                                        max = 100000, 
                                                        step = 10),
                                           numericInput("gseaParam", 
                                                        "GSEA weighting param. (rec. 1)", 
                                                        value = 1, 
                                                        min = 0, 
                                                        max = 2, 
                                                        step = 0.5),
                                           radioButtons("collapseGSEA", NULL,
                                                        choices=c("All Pathways"= "noCollapse",
                                                                  "Collapse Pathways" = "collapse")),
                                           div(id = "eaorasectui",
                                               tags$hr(),
                                               h3('ORA - selecting genes'), 
                                               numericInput("ea_subsetfdr", "P-adj cut-off", value = 0.05, min = 0, max = 1, step = 0.05),
                                               numericInput("ea_subsetlc", "LogFC cut-off", value = 1, min = 0, max = 3, step = 1),
                                               selectInput("ea_subsettype", 
                                                           "Gene status", 
                                                           c("Upregulated",
                                                             "Downregulated"), 
                                                           multiple = FALSE)
                                           ),
                                           div(id = "eagsearankingui",
                                               tags$hr(),
                                               h3('Ranking Method'), 
                                               selectInput("earankingmethod", 
                                                           "Select the ranking method (Use default for GSEA)", 
                                                           c("-log10(P-value) * sig(log2FC)",
                                                             "log Fold Change",
                                                             "-log10(P-value) * log2FC"), 
                                                           multiple = FALSE)
                                           ),
                                           
                                           
                                           actionButton("enrichementbt", "Perform analysis", class = "btn-primary")
                                  ),
                                  dropdown(label = "Plot options",
                                           icon = icon("image"),
                                           size =  "sm",
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           
                                           selectInput("ea_plottype", 
                                                       "Plot type", 
                                                       c("Dot plot",
                                                         "GSEA table plot",
                                                         "Running score and preranked list",
                                                         "Ranked list of genes"), 
                                                       multiple = FALSE),
                                           numericInput("ea_nb_categories", "Number of categories", value = 10, min = 2, max = 30, step = 1),
                                           selectInput("gsea_gene_sets", "Plot gene sets", NULL, multiple = TRUE)
                                  ),
                                  dropdown(
                                    label = "Export image",
                                    icon = icon("save"),
                                    size =  "sm",
                                    style = "bordered", 
                                    status = "primary", 
                                    width = "300px",
                                    animate = animateOptions(
                                      enter = animations$fading_entrances$fadeInLeftBig,
                                      exit = animations$fading_exits$fadeOutLeft
                                    ),
                                    tooltip = tooltipOptions(title = "Export image"),
                                    textInput("enrichementPlot.filename", label = "Filename", value = "enrichement_plot.pdf"),
                                    bsTooltip("enrichementPlot.filename", "Filename (pdf, png, svg)", "left"),
                                    numericInput("ea_width", "Figure width (in)", value = 10, min = 5, max = 30, step = 1),
                                    numericInput("ea_height", "Figure height (in)", value = 10, min = 5, max = 30, step = 1),
                                    downloadButton('saveenrichementpicture', 'Export figure',class = "btn-primary")
                                  ),
                                  dropdown(
                                    label = "Generate report",
                                    size =  "sm",
                                    icon = icon("file-code"),
                                    style = "bordered", 
                                    status = "primary", 
                                    width = "300px",
                                    animate = animateOptions(
                                      enter = animations$fading_entrances$fadeInLeftBig,
                                      exit = animations$fading_exits$fadeOutLeft
                                    ),
                                    downloadButton('reportEA', 'Download HTML report',class = "btn-primary"))
                                  ,
                                  dropdown(
                                    label = "Help material",
                                    icon = icon("question-circle"),
                                    style = "bordered", 
                                    size =  "sm",
                                    status = "primary", 
                                    width = "300px",
                                    animate = animateOptions(
                                      enter = animations$fading_entrances$fadeInLeftBig,
                                      exit = animations$fading_exits$fadeOutLeft
                                    ),
                                    shiny::actionButton(inputId='ab1', 
                                                        label = "Learn More", 
                                                        icon = icon("th"), 
                                                        onclick = "window.open('https://guangchuangyu.github.io/pathway-analysis-workshop/', '_blank')", 
                                                        class = "btn-primary"),
                                    shiny::actionButton(inputId = 'ab1', 
                                                        label="MSigDB Collections", 
                                                        icon = icon("th"), 
                                                        onclick = "window.open('http://software.broadinstitute.org/gsea/msigdb/collection_details.jsp', '_blank')", 
                                                        class = "btn-primary"))
                           ),
                           column(8,
                                  
                                  tabsetPanel(type = "pills",
                                              tabPanel("Plots",
                                                       jqui_resizable(
                                                         plotOutput("plotenrichment", height = "600")
                                                       ) #%>% withSpinner(type = 6)
                                              ),
                                              tabPanel("Table",
                                                       DT::dataTableOutput('tbl.analysis') # %>% withSpinner(type = 6)
                                              )
                                  )
                           )
                  ),
                  tabPanel("DEG Comparison Between Treatments", 
                           icon = icon("image"),
                           column(2,
                                  #sidebarPanel(id="sidebar",width = 3,
                                  dropdown(label = "DEG Output Upload",
                                           icon = icon("file-csv"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeft,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           # Input: Select a file ----
                                           fileInput("DEfile1", "Choose CSV File 1",
                                                     multiple = TRUE,
                                                     accept = c("text/csv",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv")),
                                           textInput("DEfile1Name", "Name for venn 1"),
                                           fileInput("DEfile2", "Choose CSV File 2",
                                                     multiple = TRUE,
                                                     accept = c("text/csv",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv")),
                                           textInput("DEfile2Name", "Name for venn 2"),
                                           fileInput("DEfile3", "Choose CSV File 3",
                                                     multiple = TRUE,
                                                     accept = c("text/csv",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv")),
                                           textInput("DEfile3Name", "Name for venn 3"),
                                           fileInput("DEfile4", "Choose CSV File 4",
                                                     multiple = TRUE,
                                                     accept = c("text/csv",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv")),
                                           textInput("DEfile4Name", "Name for venn 4"),
                                           fileInput("DEfile5", "Choose CSV File 5",
                                                     multiple = TRUE,
                                                     accept = c("text/csv",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv")),
                                           textInput("DEfile5Name", "Name for venn 5"),
                                           tags$div(
                                             HTML(paste(help_text3))
                                           )
                                  ),
                                  dropdown(label = "Gene Selection Criteria",
                                           icon = icon("filter"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           numericInput("vennLog2FoldChange", "log2FoldChange  cut-off:", value = 2, min = 0, max = 10, step = 0.1),
                                           numericInput("vennPadj", "P adjusted cut-off or s-value:", 0.01, min = 0, max = 1,step = 0.1)
                                  ),
                                  dropdown(label = "Venn Diagram Options",
                                           icon = icon("chart-pie"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           actionButton("venn", "Make Venn Diagram", class = "btn-primary")
                                  )
                           ),
                           column(10,
                                  tabsetPanel(type = "pills",
                                              id = "DEG",
                                              tabPanel("Venn Results",
                                                       tags$hr(),
                                                       DT::dataTableOutput('venn.results') 
                                              ),
                                              tabPanel("Venn Diagram",
                                                       tags$hr(),
                                                       plotlyOutput('vennplot') %>% withSpinner(type = 6)
                                              )
                                  )
                           )
                  ),
                  #Need to make your own custom help documents
                  tabPanel("KEGG Mapping", 
                           icon = icon("brain"),
                           column(2,
                                  #sidebarPanel(id="sidebar",width = 3,
                                  dropdown(label = "DEA results upload      ",
                                           icon = icon("file-csv"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           downloadButton('test', 'Download example DEA file',class = "btn-primary"),
                                           fileInput("testtest", "Choose CSV File",
                                                     multiple = TRUE,
                                                     accept = c("text/csv",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv"))
                                  ),
                                  
                                  #Need to add different options for the CLUEGO, pathway enrichment, settings,
                                  #also need to be able to upload custom GMT files and select necessary enrichment
                                  #files for input into CLUEGO
                                  dropdown(label = "KEGG Mapping",
                                           icon = icon("calculator"),
                                           style = "bordered", 
                                           size =  "sm",
                                           status = "primary", 
                                           width = "300px",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           selectInput("keggOrg", 
                                                       "Select your organism from KEGG", 
                                                       c("Chlamydomonas reinhardtii" = "chlam",
                                                         "Neurospora crassa" = "neur",
                                                         "Aspergillus nidulans" = "asp"), 
                                                       multiple = FALSE),
                                           numericInput("KEGGfdr", 
                                                        "P-value cut-off:", 
                                                        value = 0.05, 
                                                        min = 0, 
                                                        max = 1, 
                                                        step = 0.05),
                                           numericInput("KEGGlogFC", 
                                                        "Log2FC Cut-off:", 
                                                        value = 2, 
                                                        step = 0.5)
                                  )
                           ),
                           column(8,
                                  
                                  tabsetPanel(type = "pills",
                                              tabPanel("Plots",
                                                       jqui_resizable(
                                                         imageOutput("plotKEGG", height = "600")
                                                       ) #%>% withSpinner(type = 6)
                                              )

                                  )
                           )

                  ),

                  #Need to make your own custom help documents
                  tabPanel("Help", ##changing from tab 2, but still usibg tab2 in other parts of code
                           #icon = icon("object-group"),
                           icon = icon("question-circle"),
                           tabsetPanel(type = "pills",
                                       tabPanel("Vignette", 
                                                icon = icon("book"),
                                                includeMarkdown("Genavi.Rmd")
                                       ),
                                       tabPanel("Tutorial",
                                                icon = icon("book"),
                                                includeHTML("GENAVi_Tutorial.html")
                                       ),
                                       tabPanel("References",
                                                icon = icon("book"),
                                                includeMarkdown("References.Rmd")
                                       )
                           )
                  )
                )
)