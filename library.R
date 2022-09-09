list.of.packages <- c("Glimma",
                      "shiny",
                      "regionReport",
                      "tidyr",
                      "DESeq2",
                      "dplyr",
                      "readxl",
                      "openxlsx",
                      "iheatmapr",
                      "edgeR",
                      "ggplot2",
                      "readr",
                      "shinythemes",
                      "shinyBS",
                      "magrittr",
                      "shinyjs",
                      "shinycssloaders",
                      "shinyWidgets",
                      "prettydoc",
                      "DOSE",
                      "shinyjqui",
                      "shinydashboard",
                      "apeglm",
                      "clusterProfiler",
                      "DT",
                      "ashr",
                      "enrichplot",
                      "plotly",
                      "tidyverse",
                      "TCC",
                      "ellipse",
                      "factoextra",
                      "DEFormats",
                      "gprofiler2",
                      "fgsea",
                      "hrbrthemes",
                      "tm",
                      "proustr",
                      "VennDiagram",
                      "ICS",
                      "RColorBrewer",
                      "rrcov")

#checking missing packages from list
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

#install missing ones
if(length(new.packages)) {
  install.packages(new.packages, dependencies = TRUE,quiet = TRUE)
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    library(BiocManager)
    options(repos = BiocManager::repositories())
    BiocManager::install(new.packages, dependencies = TRUE,quiet = TRUE)
  }
  
}

suppressPackageStartupMessages({
  library(Glimma)
  library(TCC)
  library(shiny)
  library(regionReport)
  library(tidyr)
  library(DESeq2)
  library(dplyr)
  library(readxl)
  library(openxlsx)
  library(iheatmapr)
  library(edgeR)
  library(ggplot2)
  library(readr)
  library(shinythemes)
  library(shinyBS)
  library(magrittr)
  library(shinyjs)
  library(shinycssloaders)
  library(shinyWidgets)
  library(prettydoc)
  library(DOSE)
  library(shinyjqui)
  library(shinydashboard)
  library(apeglm)
  library(clusterProfiler)
  library(DT)
  library(ashr)
  library(enrichplot)
  library(plotly)
  library(tidyverse)
  library(ellipse)
  library(factoextra)
  library(DEFormats)
  library(gprofiler2)
  library(fgsea)
  library(hrbrthemes)
  library(tm)
  library(proustr)
  library(VennDiagram)
  library(RColorBrewer)
  library(ICS)
  library(rrcov)
})

# default input file
all_cell_lines <- read_csv("Chlam_Counts.csv", col_types = readr::cols(), col_names = TRUE) ##read in data going to view, try changing to all_cell_lines_ordered.csv to see if it breaks everything...no changes...change this and GENAVi.R code to change rda file

# Load gene information data
#NEED to remove this!! Not useful for our data!!
load("hg38.rda") # human
load("GRCm38.rda") # mouse

transforms <- c("raw counts", 
                "row normalized", 
                "logCPM - log Counts per Million", 
                "vst - Variance Stabilizing Transformation",
                "rlog - regularized logarithm") 

sortby <- c("-no selection-",
            "mean", 
            "standard deviation")
cell.line.clusters <- c("All genes", 
                        "Selected genes") ##do this for cell line cluster heatmaps, changed all to filtered, should "all" be an option?

custom_pal_blues <- c("#FFFFFF","#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5","#08519C","#08306B") ##try including this here in aux to use in app...

# This will be used to parse the text areas input
# possibilities of separation , ; \n
parse.textarea.input <- function(text){
  sep <- NULL
  if(grepl(";",text)) sep <- ";"
  if(grepl(",",text)) sep <- ","
  if(grepl("\n",text)) sep <- "\n"
  if(is.null(sep)) {
    text <- text
  } else {
    text <- unlist(stringr::str_split(text,sep))
  }
  return (text)
}


createTable2 <- function(df,tableType = "HOMRNA", show.rownames = TRUE){
  DT::datatable(df,
                extensions = c('Buttons',"FixedHeader"),
                class = 'cell-border stripe',
                options = list(dom = 'Blfrtip',
                               initComplete = JS(
                                 "function(settings, json) {",
                                 "$(this.api().table().header()).css({'background-color': '#F0F8FF', 'color': '#000000'});",
                                 "}"),
                               buttons =
                                 list('colvis', list(
                                   extend = 'collection',
                                   buttons = list(list(extend='csv',
                                                       filename = tableType),
                                                  list(extend='excel',
                                                       filename = tableType),
                                                  list(extend = 'pdf',
                                                       title = "",
                                                       filename = tableType)),
                                   text =  'Download current view'
                                 )),
                               fixedHeader = FALSE,
                               pageLength = 20,
                               scrollX = TRUE,
                               lengthMenu = list(c(10, 20, 50,100), c('10', '20', '50','100'))
                ),
                rownames = show.rownames,
                filter   = 'top'
  )
}

createTable <- function(df,selected_rows=NULL,tableType = "HOMRNA", show.rownames = FALSE){
  DT::datatable(df,
                extensions = c('Buttons',"FixedHeader"),
                class = 'cell-border stripe',
                selection = list(mode = "multiple", target= 'row', selected = selected_rows),
                options = list(dom = 'Blfrtip',
                               initComplete = JS(
                                 "function(settings, json) {",
                                 "$(this.api().table().header()).css({'background-color': '#F0F8FF', 'color': '#000000'});",
                                 "$(this.api().table().header()).css({'color': '#000000'});",
                                 "}"),
                               columnDefs = list(
                                 list(visible=FALSE, targets=c(0))
                               ),
                               order = c(0,"desc"),
                               deferRender = TRUE,
                               paging = T,
                               buttons =
                                 list('colvis', list(
                                   extend = 'collection',
                                   buttons = list(list(extend='csv',
                                                       filename = tableType),
                                                  list(extend='excel',
                                                       filename = tableType),
                                                  list(extend='pdf',
                                                       title = "",
                                                       filename= tableType)),
                                   text = 'Download current view'
                                 )),
                               fixedHeader = FALSE,
                               pageLength = 20,
                               scrollX = TRUE,
                               lengthMenu = list(c(10, 20, 50,100), c('10', '20', '50','100'))
                ),
                rownames = show.rownames,
                filter   = 'top'
  )
}

edgerConts <- function(levels, reference) {
  vector <- c()
  lev <- levels[levels != reference]
  for (i in lev) {
    vector[i] <- paste(i,reference, sep="-")
  }
  return(vector)
}

rownorm <- function(counts.filtered)
{
  rownorm.tbl <-  (counts.filtered - rowMeans(counts.filtered,na.rm = TRUE)) / apply(counts.filtered,1,sd)
  colnames(rownorm.tbl) <- colnames(counts.filtered)
  rownorm.tbl
}

getEndGeneInfo <- function(data){
  
  numeric.cols <-   sum(sapply(colnames(data), function(x) {unique(class(data[[x]]) %in% c("integer","numeric"))}))
  if(ncol(data) - numeric.cols > 2) { # we have metadata
    idx <- sum(ifelse(grepl("Start",colnames(data),ignore.case = T),class(data[[grep("Start",colnames(data),ignore.case = T)]])  %in% c("numeric","integer"),0),
               ifelse(grepl("End",colnames(data),ignore.case = T),class(data[[grep("End",colnames(data),ignore.case = T)]]) %in% c("numeric","integer"),0),
               ifelse(grepl("Chr",colnames(data),ignore.case = T),class(data[[grep("Chr",colnames(data),ignore.case = T)]]) %in% c("numeric","integer"),0),
               ifelse(grepl("Strand",colnames(data),ignore.case = T),class(data[[grep("Strand",colnames(data),ignore.case = T)]])  %in% c("numeric","integer"),0),
               ifelse(grepl("Length",colnames(data),ignore.case = T),class(data[[grep("Length",colnames(data),ignore.case = T)]])  %in% c("numeric","integer"),0))
    ngene <- ncol(data) - numeric.cols + idx
  } else {
    nsamples <- ncol(data) - 1 
    data <- addgeneinfo(data) 
    ngene <-  ncol(data) - nsamples 
  }
  return(list(data = data,ngene = ngene))
}

addgeneinfo <- function(data){
  id <- data %>% pull(1)
  
  mouse.only.symbols <- setdiff(GRCm38$Symbol,hg38$Symbol)
  human.only.symbols <- setdiff(hg38$Symbol,GRCm38$Symbol)
  
  if(all(grepl("ENSG",id))) {
    colnames(data)[1] <- "EnsemblID"
    aux <- strsplit(data$EnsemblID,"\\.")
    data$EnsemblID <- as.character(unlist(lapply(aux,function(x) x[1])))
    data <- merge(hg38,data,by = "EnsemblID")
  } else  if(all(grepl("ENSMUSG",id))) { # Mouse
    colnames(data)[1] <- "EnsemblID"
    aux <- strsplit(data$EnsemblID,"\\.")
    data$EnsemblID <- as.character(unlist(lapply(aux,function(x) x[1])))
    data <- merge(GRCm38,data,by = "EnsemblID")
  } else {
    colnames(data)[1] <- "Symbol"
    data <- data
  }
  return(data)
}

help_text <-  '<div class="panel panel-primary">
  <div class="panel-heading"> <span style="padding-left:10px"><b> Input file description</b> </span></div>
<div class="panel-body">
<style type="text/css">
.tg {
border-collapse: collapse;
border-spacing: 0;
border: none;
}
.tg td {
font-family: Arial, sans-serif;
font-size: 14px;
padding: 10px 5px;
border-style: solid;
border-width: 0px;
overflow: hidden;
word-break: normal;
}
.tg th {
font-family: Arial, sans-serif;
font-size: 14px;
font-weight: normal;
padding: 10px 5px;
border-style: solid;
border-width: 0px;
overflow: hidden;
word-break: normal;
}
.tg .tg-s6z2 {
text-align: center
}
</style>
<table class="tg">
<tr>
<th class="tg-031e"> <span class="label label-primary"> Format</span></th>
<th class="tg-031e"> comma-separated values (CSV)
</tr>
<tr>
<th class="tg-031e"> <span class="label label-primary"> Column 1</span></th>
<th class="tg-031e"> Gene Identifier -  Gene symbol or ENSEMBL ID
</tr>
<tr>
<th class="tg-031e"> <span class="label label-primary"> Column 2-n</span></th>
<th class="tg-031e"> Gene expression raw counts </th>
</tr>
<tr>
</table>
</div>
</div>'

help_text2 <-  '<div class="panel panel-primary">
  <div class="panel-heading"> <span style="padding-left:10px"><b> Input file description</b> </span></div>
<div class="panel-body">
<style type="text/css">
.tg {
border-collapse: collapse;
border-spacing: 0;
border: none;
}
.tg td {
font-family: Arial, sans-serif;
font-size: 14px;
padding: 10px 5px;
border-style: solid;
border-width: 0px;
overflow: hidden;
word-break: normal;
}
.tg th {
font-family: Arial, sans-serif;
font-size: 14px;
font-weight: normal;
padding: 10px 5px;
border-style: solid;
border-width: 0px;
overflow: hidden;
word-break: normal;
}
.tg .tg-s6z2 {
text-align: center
}
</style>
<table class="tg">
<tr>
<th class="tg-031e"> <span class="label label-primary"> Format</span></th>
<th class="tg-031e"> comma-separated values (CSV)
</tr>
<tr>
<th class="tg-031e"> <span class="label label-primary"> Column 1</span></th>
<th class="tg-031e"> Sample ID (same from raw counts data)
</tr>
<tr>
<th class="tg-031e"> <span class="label label-primary"> Column 2-n</span></th>
<th class="tg-031e"> Other metadata (condition, covariates) </th>
</tr>
<tr>
</table>
</div>
</div>'

help_text3 <-  '<div class="panel panel-primary">
  <div class="panel-heading"> <span style="padding-left:10px"><b> Input file description</b> </span></div>
<div class="panel-body">
<style type="text/css">
.tg {
border-collapse: collapse;
border-spacing: 0;
border: none;
}
.tg td {
font-family: Arial, sans-serif;
font-size: 14px;
padding: 10px 5px;
border-style: solid;
border-width: 0px;
overflow: hidden;
word-break: normal;
}
.tg th {
font-family: Arial, sans-serif;
font-size: 14px;
font-weight: normal;
padding: 10px 5px;
border-style: solid;
border-width: 0px;
overflow: hidden;
word-break: normal;
}
.tg .tg-s6z2 {
text-align: center
}
</style>
<table class="tg">
<tr>
<th class="tg-031e"> <span class="label label-primary"> Format</span></th>
<th class="tg-031e"> comma-separated values (CSV)
</tr>
<tr>
<th class="tg-031e"> <span class="label label-primary"> Column 1</span></th>
<th class="tg-031e"> Gene ID
<tr>
</table>
</div>
</div>'