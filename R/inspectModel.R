inspectModel <- function(model, runShiny = TRUE){
  modelClone <- clonePsydiffModel(model)
  filename <- tempfile(pattern = "modelClone", tmpdir = tempdir(), fileext = ".RData")
  save(modelClone, file = filename)
  parameters <- getParameterValues(modelClone)
  parameters <- parameters[sort(names(parameters))]

  UI <- paste0('
   library(shiny)
   library(ggplot2)
   library(gridExtra)
   library(psydiff)
   load("', filename, '")
   parameters <- getParameterValues(modelClone)

   ui <- fluidPage(
  titlePanel("Psydiff model inspection"),
  sidebarPanel(
  actionButton("close", "Close and return values"),
  checkboxInput("simulateOnly", "Simulate only", value = FALSE, width = NULL),
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
  parameterValues <- getParameterValues(modelClone)
  simulateOnly <- as.logical(input$simulateOnly)
', SERVER, '
    setParameterValues(parameterTable = modelClone$pars$parameterTable,
                       parameterValues = parameterValues,
                       parameterLabels = names(parameterValues))
    if(simulateOnly){
      simDat <- simulateData(modelClone)
      predictedManifest <- simDat$predictedManifest
      plotTitle <- "Simulated paths for observed data"
    }else{
      f <- fitModel(modelClone)
      predictedManifest <- f$predictedManifest
      plotTitle <- paste0("Predicted observations; m2LL = ", sum(f$m2LL))
    }
observedManifest <- modelClone$data$observations
colnames(predictedManifest) <- paste0("Y", seq_len(modelClone$nmanifest))
colnames(observedManifest) <- paste0("Y", seq_len(modelClone$nmanifest))

nrowsOfPlot <- floor(modelClone$nmanifest/2) +1
par(mfrow = c(nrowsOfPlot, 2))

times <- c()
for(p in unique(modelClone$data$person)){
  times <- c(times, cumsum(modelClone$data$dt[modelClone$data$person == p]))
}

df <- data.frame("person" = as.factor(modelClone$data$person), "times" = times,  predictedManifest)
obs <- data.frame("person" = as.factor(modelClone$data$person), "times" = times,  observedManifest)
')
plots <- paste0('plot', seq_len(modelClone$nmanifest), ' <- ggplot(data = df, aes(x=times, y= Y', seq_len(modelClone$nmanifest),', group=person, color=person)) +
           geom_line()+
    geom_point(data = obs,
               mapping = aes(x = times, y = Y', seq_len(modelClone$nmanifest),', color = person)) +
                ggtitle(plotTitle)', collapse = " \n")
arrangePlots <- paste0('
    gr = grid.arrange(',paste0('plot', seq_len(modelClone$nmanifest), collapse = ", "),', ncol= ',ifelse(modelClone$nmanifest>1,2,1),')
                       print(gr)')
endSERVER <- '
  });
  observe({
      if(input$close > 0){
        stopApp(getParameterValues(modelClone))
      }
  })
  }
shinyApp(ui, server)
'
combined <- paste0(UI, SERVER, plots, arrangePlots, endSERVER)

if(runShiny){
  filename <- tempfile(pattern = "psydiff_Shiny", tmpdir = tempdir(), fileext = ".R")
  fileConn<-file(filename)
  writeLines(combined, fileConn)
  close(fileConn)

  retValues <- shiny::runApp(filename)
  return(retValues)
}else{
  return(combined)
}
}

