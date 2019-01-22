jscode <- 
"
shinyjs.disableTab = function(name) {

var tab = $('.nav li a[data-value=' + name + ']');
tab.bind('click.tab', function(e) {
e.preventDefault();
return false;
});
tab.addClass('disabled');
}

shinyjs.enableTab = function(name) {
var tab = $('.nav li a[data-value=' + name + ']');
tab.unbind('click.tab');
tab.removeClass('disabled');
}
"

css <-
"
.nav li a.disabled {
background-color: #ccc !important;
color: #aaa !important;
cursor: not-allowed !important;
border-color: #aaa !important;
}
"

library(shiny)
library(shinyjs)
#library(DT)
source("iCAT.R")
options(warn=1)

ui <- fluidPage(useShinyjs(),
                extendShinyjs(text = jscode, functions =c("disableTab", "enableTab")),
                inlineCSS(css),
                
                img(src='cat2.png', height=120, width=160, hspace = 5, vspace = 5),
                
                tabsetPanel(id = 'navbar',
                  tabPanel(title = "Training",
                           id = "trnTab",
                           fluid = T,
                           sidebarLayout(
                             sidebarPanel(
                               fileInput(
                                 "pre",
                                 "Negative Training",
                                 multiple = TRUE,
                                 accept = c("text/tsv", "text/tab-seperated-values", ".tsv")
                               ),
                               fileInput(
                                 "post",
                                 "Positive Training",
                                 multiple = TRUE,
                                 accept = c("text/tsv", "text/tab-seperated-values", ".tsv")
                               ),
      
                               radioButtons(
                                 "field",
                                 "Analyze Clonotypes By:",
                                 choices = list(
                                   "CDR3 Amino Acid Only" = "aminoAcid",
                                   "TCRV-CDR3-TCRJ" = "vGeneName aminoAcid jGeneName",
                                   "Nucleic Acid (DNA)" = "nucleotide"
                                 )
                               ),
                               textInput("pcut", "Max p-value", value = 0.1),
                               textInput("thresh", "Min Threshold of Public Sequences", value = 1),
                               tags$hr(),
                               
                               actionButton("run", "Train Model"),
                               br(),
                               br(),
                               div(hidden(downloadButton('dnScreen', label="Save Parameters", style='padding-left:125px; padding-right:125px')))
                             ),
                             
                             mainPanel(
                               hidden(h4(id = 'h4', "Training Data Summary:")),
                               tableOutput('trnTable'),
                               hidden(downloadButton('dnSummary', label="Table")),
                               hidden(h4(id = 'h1', "Pre/Post Distributions:")),
                               plotOutput(outputId = "plot"),
                               div(style="display:inline-block", hidden(downloadButton('dnPlot', label="Plot PNG"))),
                               div(style="display:inline-block", hidden(downloadButton('dnPlotPDF', label="Plot PDF"))),
                               br(),
                               br(),
                               hidden(h4(id = 'h2', "Classification Matrix: ")),
                               tableOutput('table'),
                               hidden(downloadButton('dnClass', label="Table"))
                             )
                             
                           )
                           ),
                  
                  tabPanel(
                    title = "Library",
                    value = "libTab",
                    DT::dataTableOutput("library"),
                    br(),
                    hidden(downloadButton('dnLib', label="Table"))
                  ),
                  
                  tabPanel(
                    title = "Prediction",
                    value = "predTab",
                    fluid = T,
                    sidebarLayout(
                      sidebarPanel(
                        fileInput(
                          "indpt",
                          "Independent Sample(s)",
                          multiple = T,
                          accept = c("text/tsv", "text/tab-separated-values", ".tsv")
                        ),
                        tags$hr(),
                        actionButton("pred", "Predict Independent Sample(s)")
                      ),
                      mainPanel(hidden(h4(
                        id = "h3", "Prediction Results:"
                      )),
                      DT::dataTableOutput('result'),
                      br(),
                      hidden(downloadButton('dnPred', label="Table"))
                      )
                    )
                  )
                )
)

server <- function(input, output, session) {
  js$disableTab("predTab")
  js$disableTab("libTab")
  
  options(shiny.maxRequestSize=3000*1024^2)
  
  both <- eventReactive(input$run, {
    validate(need(input$pre !="", "Please select negative samples"),
               need(input$post != "", "Please select postive samples" ))
    progress <- shiny::Progress$new(style = 'notification')
    progress$set(message = "Training Procedure", value = 0)
    on.exit(progress$close())
    
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <-  value + (progress$getMax() - value) / 5
      }
      progress$set(value = value, detail = detail)
    }
    
    updateProgress(detail = "Reading files")
    naive <- readPre(input$pre$datapath, input$field)
    vaccs <- readPost(input$post$datapath, input$field)
    
    anl <- analyse(naive, vaccs, input$pre$datapath, input$post$datapath, input$field, input$pcut, input$thresh, updateProgress)
    show("h1")
    show("h2")
    show("h4")
    show("dnPlot")
    show("dnPlotPDF")
    show("dnLib")
    show("dnScreen")
    show("libTab")
    show("dnSummary")
    show("dnClass")
    
    js$enableTab("predTab")
    js$enableTab("libTab") 
    return(anl)
  })
  
  preds <- eventReactive(input$pred, {
    show("h3")
    show("dnPred")
    return(pred(both(), readPost(input$indpt$datapath, input$field), input$indpt$datapath, input$indpt$name, input$field))
  })
  
  output$plot <- renderPlot({
    plotHist(both())
  })
  
  output$trnTable <- renderTable({
    both()
    trnStats(input$pre$datapath, input$post$datapath, input$field)
  }, rownames = T)
  output$table <- renderTable({
    classMat(both())
  },
  include.rownames = T)
  output$result  <- DT::renderDataTable({
    DT::datatable(preds()) %>%
      DT::formatStyle(
        'Prediction',
        target = 'row',
        color = DT::styleEqual(c("Negative", "Positive"), c('blue', 'red'))
      )
  })
  
  output$dnPred <- downloadHandler(
    filename = "predictions.csv",
    content = function(file) {
      write.csv(preds(), file, row.names=F)
    }
  )
  output$dnClass <- downloadHandler(
    filename = "classification_matrix.csv",
    content = function(file) {
      write.csv(classMat(both()), file, row.names=T)
    }
  )
  output$dnSummary <- downloadHandler(
    filename = "training_summary.csv",
    content = function(file) {
      both()
      write.csv(trnStats(input$pre$datapath, input$post$datapath, input$field), file, row.names=T)
    }
  )
  
  output$dnPlot <- downloadHandler(
    filename = "clonotypes_distribution.png",
    content = function(file) {
      both()
      png(file)
      print(plotHist(both()))
      dev.off()
    }
  )
  output$dnScreen <- downloadHandler(
    filename = paste("iCAT_report_", Sys.time(), ".txt", sep = ""),
    content = function(file) {
      cat(paste("iCAT Run on ", Sys.time()), file=file, sep="\n")
      cat("\nNegative Files:", file=file, append=T, sep="\n")
      cat(input$pre$name, file=file, append=T, sep="\n")
      cat("\nPositive Files:", file=file, append=T, sep="\n")
      cat(input$post$name, file=file, append=T, sep="\n")
      cat("\nField:", file=file, append=T, sep="\n")
      cat(input$field, file=file, append=T, sep="\n")
      cat("\nMax PValue:", file=file, append=T, sep="\n")
      cat(input$pcut, file=file, append=T, sep="\n")
      cat("\nMin Threshold of Public Sequences:", file=file, append=T, sep="\n")
      cat(input$thresh, file=file, append=T, sep="\n")
      
    }
  )
  output$dnPlotPDF <- downloadHandler(
    filename = "clonotypes_distribution.pdf",
    content = function(file) {
      pdf(file)
      print(plotHist(both()))
      dev.off()
    }
  )
  output$dnLib <- downloadHandler(
    filename = "clonotypes_library.csv",
    content = function(file) {
      write.csv(getLib(), file, row.names=F)
    }
  )

  output$library <- DT::renderDataTable({
    both()
    DT::datatable(getLib())
  },
  
  options = list(scrollX = TRUE))
}

shinyApp(ui, server)


