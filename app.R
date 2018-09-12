library(shiny)
source("iCAT.R")
options(warn=1)

ui <- fluidPage(
  titlePanel("iCAT Web App"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("pre",
                "Choose PRE directory:",
                value = paste(getwd(), "/Pre/", sep = "")
      ),
      textInput("post",
                "Choose POST directory:",
                value = paste(getwd(), "/Post/", sep = "")
      ),
      tags$hr(),
      
      actionButton("run", "Run Analysis")
    ),
    
    mainPanel(
      plotOutput(outputId = "plot")
    )
    
  )
)

server <- function(input, output) {
  
  both <- eventReactive(input$run, {
    lpre <- listDir(input$pre)
    lpost <- listDir(input$post)
    naive <- readPre(lpre)
    vaccs <- readPost(lpost)
    t1 <- Sys.time()
    anl <- analyse(naive, vaccs, lpre, lpost)
    t2 <- Sys.time()
    print(round(difftime(t2, t1, units = "secs"), digits = 2))
    return(anl)
  })
  
  output$plot <- renderPlot({
    plotHist(both())
  })
}

shinyApp(ui, server)
