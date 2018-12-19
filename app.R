library(shiny)
#library(DT)
source("iCAT.R")
options(warn=1)

ui <- fluidPage(
  titlePanel("iCAT Web App"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("pre",
                "Negative Training", multiple = TRUE,
                accept = c("text/tsv", "text/tab-seperated-values", ".tsv")
      ),
      fileInput("post",
                "Positive Training", multiple = TRUE,
                accept = c("text/tsv", "text/tab-seperated-values", ".tsv")
      ),
      fileInput("indpt",
                "Independent", multiple =T,
                accept = c("text/tsv", "text/tab-seperated-values", ".tsv")),
      radioButtons("field", "Analyze Clonotypes By:",
                   choices= list("CDR3 Amino Acid Only" = "aminoAcid",
                                 "TCRV-CDR3-TCRJ" = "aminoAcid vGeneName jGeneName",
                                 "Nucleic Acid (DNA)" = "nucleotide")),
      textInput("pcut", "Max p-value", value = 0.1),
      textInput("thresh", "Min Threshold of Public Sequences", value = 1),
      tags$hr(),
      
      actionButton("run", "Train Model"),
      actionButton("pred", "Predict Independent Sample")
    ),
    
    mainPanel(
      tags$h3("Pre/Post Distributions: "),
      plotOutput(outputId = "plot"),
      br(),
      br(),
      tags$h3("Classification Matrix: "),
      tableOutput('table'),
      br(),
      br(),
      tags$h3("Independent Prediction: "),
      textOutput('result')
    )
    
  )
)

server <- function(input, output) {
  options(shiny.maxRequestSize=30*1024^2)
  
  both <- eventReactive(input$run, {
    naive <- readPre(input$pre$datapath, input$field)
    vaccs <- readPost(input$post$datapath, input$field)
    t1 <- Sys.time()
    anl <- analyse(naive, vaccs, input$pre$datapath, input$post$datapath, input$field, input$pcut, input$thresh)
    
    t2 <- Sys.time()
    print(round(difftime(t2, t1, units = "secs"), digits = 2))
    return(anl)
  })
  
  preds <- eventReactive(input$pred, {
    
    return(pred(both(), readPost(input$indpt$datapath, input$field), input$indpt$datapath, input$field))
  })
  
  output$plot <- renderPlot({
    plotHist(both())
  })
  
  output$table <- renderTable({
    classMat(both())
  })
  output$result  <- renderText({
    preds()
  })
}

shinyApp(ui, server)


