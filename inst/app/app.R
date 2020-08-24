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
ui <- shiny::fluidPage(shinyjs::useShinyjs(),
                shinyjs::extendShinyjs(text = jscode, functions =c("disableTab", "enableTab")),
                shinyjs::inlineCSS(css),

                shiny::img(src='cat2.png', height=120, width=160, hspace = 5, vspace = 5),

                shiny::tabsetPanel(id = 'navbar',
                  shiny::tabPanel(title = "Training",
                           id = "trnTab",
                           fluid = TRUE,
                           shiny::sidebarLayout(
                             shiny::sidebarPanel(
                               shiny::fileInput(
                                 "pre",
                                 "Negative Training",
                                 multiple = TRUE,
                                 accept = c("text/tsv", "text/tab-seperated-values", ".tsv")
                               ),
                               shiny::fileInput(
                                 "post",
                                 "Positive Training",
                                 multiple = TRUE,
                                 accept = c("text/tsv", "text/tab-seperated-values", ".tsv")
                               ),
                               
                               shiny::radioButtons(
                                 "field",
                                 "Analyze Clonotypes By (column names in parantheses):",
                                 choices = list(
                                   "CDR3 Amino Acid (aminoAcid)" = "aminoAcid",
                                   "TCRV-CDR3-TCRJ (vGeneName aminoAcid jGeneName)" = "vGeneName aminoAcid jGeneName",
                                   "Nucleic Acid (nucleotide)" = "nucleotide",
                                   "Other" = "other"
                                 )
                               ),
                               shinyjs::hidden(shiny::textInput("otherbt", "Please provide data columns to analyze by (space-separated)",
                                                value = "vGeneName aminoAcid jGeneName")),
                               
                               shiny::radioButtons(
                                 "count",
                                 "What column name contains clonotype counts in your data?",
                                 choices = list(
                                   '"count (templates)"' = "count (templates)",
                                   '"copy"' = "copy",
                                   "Other" = "other"
                                 )
                               ),
                               shinyjs::hidden(shiny::textInput("otherCount", "Please provide the column name for clontype counts",
                                                                value = "copy")),
                               shiny::textInput("pcut", "Max p-value", value = 0.1),
                               shiny::textInput("copyrange", "Range of Acceptable Copies per Clonotype", value = "1 99"),
                               shiny::textInput("thresh", "Min # of Public Sequences", value = 1),
                               shiny::tags$hr(),

                               shiny::actionButton("run", "Train Model"),
                               shiny::br(),
                               shiny::br(),
                               shiny::div(shinyjs::hidden(shiny::downloadButton('dnScreen', label="Save Parameters", style='padding-left:125px; padding-right:125px')))
                             ),

                             shiny::mainPanel(
                               shinyjs::hidden(shiny::h4(id = 'h4', "Training Data Summary:")),
                               shiny::tableOutput('trnTable'),
                               shinyjs::hidden(shiny::downloadButton('dnSummary', label="Table")),
                               shinyjs::hidden(shiny::h4(id = 'h1', "Pre/Post Distributions:")),
                               shiny::plotOutput(outputId = "plot"),
                               shiny::div(style="display:inline-block", shinyjs::hidden(shiny::downloadButton('dnPlot', label="Plot PNG"))),
                               shiny::div(style="display:inline-block", shinyjs::hidden(shiny::downloadButton('dnPlotPDF', label="Plot PDF"))),
                               shiny::br(),
                               shiny::br(),
                               shinyjs::hidden(shiny::h4(id = 'h2', "Classification Matrix: ")),
                               shiny::tableOutput('table'),
                               shinyjs::hidden(shiny::downloadButton('dnClass', label="Table"))
                             )

                           )
                           ),

                  shiny::tabPanel(
                    title = "Library",
                    value = "libTab",
                    DT::dataTableOutput("library"),
                    shiny::br(),
                    shinyjs::hidden(shiny::downloadButton('dnLib', label="Table"))
                  ),

                  shiny::tabPanel(
                    title = "Prediction",
                    value = "predTab",
                    fluid = TRUE,
                    shiny::sidebarLayout(
                      shiny::sidebarPanel(
                        shiny::fileInput(
                          "indpt",
                          "Independent Sample(s)",
                          multiple = TRUE,
                          accept = c("text/tsv", "text/tab-separated-values", ".tsv")
                        ),
                        shiny::tags$hr(),
                        shiny::actionButton("pred", "Predict Independent Sample(s)")
                      ),
                      shiny::mainPanel(shinyjs::hidden(shiny::h4(
                        id = "h3", "Prediction Results:"
                      )),
                      DT::dataTableOutput('result'),
                      shiny::br(),
                      shinyjs::hidden(shiny::downloadButton('dnPred', label="Table"))
                      )
                    )
                  )
                )
)

server <- function(input, output, session) {
  shinyjs::js$disableTab("predTab")
  shinyjs::js$disableTab("libTab")

  options(shiny.maxRequestSize=10000*1024^2)

  observe({
    shinyjs::toggle("otherbt", anim = TRUE, condition = input$field == "other")
    if (input$field == "other") {
      field <<- input$otherbt
    } else {
      field <<- input$field
    }
  })
  
  observe({
    shinyjs::toggle("otherCount", anim = TRUE, condition = input$count == "other")
    if (input$count == "other") {
      count <<- input$otherCount
    } else {
      count <<- input$count
    }
  })

  both <- shiny::eventReactive(input$run, {
    shiny::validate(shiny::need(input$pre !="", "Please select negative samples"),
                    shiny::need(input$post != "", "Please select postive samples" ))
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
    copyrange <- input$copyrange

    updateProgress(detail = "Reading files")
    naive <- readTrn(input$pre$datapath, field, count, copyrange, "naive")
    vaccs <- readTrn(input$post$datapath, field, count, copyrange, "vacc")
    anl <- train(naive, vaccs, input$pre$datapath, input$post$datapath, field, count, copyrange, input$pcut, input$thresh, updateProgress)
    shinyjs::show("h1")
    shinyjs::show("h2")
    shinyjs::show("h4")
    shinyjs::show("dnPlot")
    shinyjs::show("dnPlotPDF")
    shinyjs::show("dnLib")
    shinyjs::show("dnScreen")
    shinyjs::show("libTab")
    shinyjs::show("dnSummary")
    shinyjs::show("dnClass")

    shinyjs::js$enableTab("predTab")
    shinyjs::js$enableTab("libTab")
    return(anl)
  })

  preds <- shiny::eventReactive(input$pred, {
    shinyjs::show("h3")
    shinyjs::show("dnPred")
    return(pred(both(), input$indpt$datapath, input$indpt$name, field, count, copyrange))
  })

  output$plot <- shiny::renderPlot({
    plotHist(both())
  })

  output$trnTable <- shiny::renderTable({
    both()
    trnStats(input$pre$datapath, input$post$datapath, field, count, input$copyrange)
  }, rownames = TRUE)
  output$table <- renderTable({
    classMat(both())
  },
  include.rownames = TRUE)
  output$result  <- DT::renderDataTable({
      DT::formatStyle(
        DT::datatable(preds()),
        'Prediction',
        target = 'row',
        color = DT::styleEqual(c("Negative", "Positive"), c('blue', 'red'))
      )
  })

  output$dnPred <- shiny::downloadHandler(
    filename = "predictions.csv",
    content = function(file) {
      write.csv(preds(), file, row.names=FALSE)
    }
  )
  output$dnClass <- shiny::downloadHandler(
    filename = "classification_matrix.csv",
    content = function(file) {
      write.csv(classMat(both()), file, row.names=TRUE)
    }
  )
  output$dnSummary <- shiny::downloadHandler(
    filename = "training_summary.csv",
    content = function(file) {
      both()
      write.csv(trnStats(input$pre$datapath, input$post$datapath, field, count, input$copyrange), file, row.names=TRUE)
    }
  )

  output$dnPlot <- shiny::downloadHandler(
    filename = "clonotypes_distribution.png",
    content = function(file) {
      both()
      png(file)
      print(plotHist(both()))
      dev.off()
    }
  )
  output$dnScreen <- shiny::downloadHandler(
    filename = paste("iCAT_report_", Sys.time(), ".txt", sep = ""),
    content = function(file) {
      cat(paste("iCAT Run on ", Sys.time()), file=file, sep="\n")
      cat("\nNegative Files:", file=file, append=TRUE, sep="\n")
      cat(input$pre$name, file=file, append=TRUE, sep="\n")
      cat("\nPositive Files:", file=file, append=TRUE, sep="\n")
      cat(input$post$name, file=file, append=TRUE, sep="\n")
      cat("\nField:", file=file, append=TRUE, sep="\n")
      cat(field, file=file, append=TRUE, sep="\n")
      cat("\nCount column name:", file=file, append=TRUE, sep="\n")
      cat(count, file=file, append=TRUE, sep="\n")
      cat("\nMax PValue:", file=file, append=TRUE, sep="\n")
      cat(input$pcut, file=file, append=TRUE, sep="\n")
      cat("\nMin Threshold of Public Sequences:", file=file, append=TRUE, sep="\n")
      cat(input$thresh, file=file, append=TRUE, sep="\n")

    }
  )
  output$dnPlotPDF <- shiny::downloadHandler(
    filename = "clonotypes_distribution.pdf",
    content = function(file) {
      pdf(file)
      print(plotHist(both()))
      dev.off()
    }
  )
  output$dnLib <- shiny::downloadHandler(
    filename = "clonotypes_library.csv",
    content = function(file) {
      write.csv(getLib(both()), file, row.names=FALSE)
    }
  )

  output$library <- DT::renderDataTable({
    DT::datatable(getLib(both()))
  },

  options = list(scrollX = TRUE))
}

shinyApp(ui, server)
