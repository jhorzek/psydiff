inspectModel <- function(model){

  StartSettings <- model
  filename <- tempfile(pattern = "StartSettings", tmpdir = tempdir(), fileext = ".RData")
  save(StartSettings, file = filename)
  parameters <- getParameterValues(model)

  UI <- paste0('
   library(shiny)
   library(psydiff)
   load("', filename, '")
   parameters <- getParameterValues(model)

   ui <- fluidPage(
  titlePanel("Psydiff model inspection"),
  sidebarPanel(
  ')

  SERVERHEAD <- '
  server <- function(input, output, session) {
  '
  SERVER <- ''
  numpar <- length(parameters)
  for(parameter in names(parameters)){
    addelemUi <- paste0('
      numericInput(
        inputId = "',parameter, '", label = "',parameter, '",
        value = ', parameters[parameter], ',
        step =  .1
        )'
    )

    UI <- paste0(UI, addelemUi, sep =  ifelse(parameter == names(parameters)[numpar], "", ","))

    addelemServer <- paste0('
      parameterValues["',parameter, '"]  = as.numeric(input$',parameter, ')
      '
    )
    SERVER <-  paste0(SERVER, addelemServer)
  }


  UI <- paste0(UI, '),
  mainPanel(
    plotOutput("simPaths")
  )
)
')


SERVER <- paste0(SERVERHEAD,'

output$simPaths <- renderPlot({
  parameterValues <- getParameterValues(model)
', SERVER, '
    setParameterValues(parameterTable = model$pars$parameterTable,
                       parameterValues = parameterValues,
                       parameterLabels = names(parameterValues))
    simDat <- simulateData(model)
    matplot(simDat$predictedManifest, type = "l")
    points(simDat$simulatedObservation)
  })
  }
shinyApp(ui, server)
')
combined <- paste0(UI, SERVER)

filename <- tempfile(pattern = "psydiff_Shiny", tmpdir = tempdir(), fileext = ".R")
fileConn<-file(filename)
writeLines(combined, fileConn)
close(fileConn)

runApp(filename)
return(combined)
}

