library(shiny)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput("counts", "Input Counts csv File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      tags$hr(),
      checkboxInput("header1", "Header", TRUE),

      fileInput("meta", "Input Metadata csv File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      tags$hr(),
      checkboxInput("header2", "Header", TRUE),
      
      actionButton(inputId= "checkHead", 
                   label="Check that columns match rows")
      
    ),
    mainPanel(
      verbatimTextOutput("checkData"),
      verbatimTextOutput("user_opt"),
      tableOutput("countcontents"),
      tableOutput("metacontents")
    )
  )
)

server <- function(input, output) {
  statement <- reactiveValues(data="Need to Check Rows and Columns Match")
  
  countData <- reactive({
    req(input$counts, input$header1, file.exists(input$counts$datapath))
    read.csv(input$counts$datapath, header = input$header1)
  })
  output$countcontents <- renderTable({
    req(head(countData()))
    head(countData())
  })
  

  metadata <- reactive({
    req(input$meta, input$header2, file.exists(input$counts$datapath))
    read.csv(input$meta$datapath, header = input$header2)
  })
  output$metacontents <- renderTable({
    req(metadata())
    metadata()
  })
  
  observeEvent(input$checkHead, {
    statement$data <- print(ncol(metadata))

  })
  
  output$checkData <- renderText({
    statement$data
  })
  
}

shinyApp(ui, server)