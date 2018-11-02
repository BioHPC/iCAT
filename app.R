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
                                 "TCRV-CDR3-TCRJ" = "b",
                                 "Nucleic Acid (DNA)" = "nucleotide")),
      textInput("pvalue", "Max p-value", value = 0.2),
      textInput("thresh", "Min Threshold of Public Sequences", value = 3),
      tags$hr(),
      
      actionButton("run", "Run Analysis")
    ),
    
    mainPanel(
      plotOutput(outputId = "plot")
    )
    
  )
)

server <- function(input, output) {
  options(shiny.maxRequestSize=30*1024^2)
  
  both <- eventReactive(input$run, {
    naive <- readPre(input$pre$datapath, input$field)
    vaccs <- readPost(input$post$datapath, input$field)
    t1 <- Sys.time()
    anl <- analyse(naive, vaccs, input$pre$datapath, input$post$datapath, input$field)
    t2 <- Sys.time()
    print(round(difftime(t2, t1, units = "secs"), digits = 2))
    return(anl)
  })
  
  output$plot <- renderPlot({
    plotHist(both())
  })
}

shinyApp(ui, server)
