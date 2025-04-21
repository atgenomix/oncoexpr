#' @title PCA plot UI
#' @description PCA
#' @param id module id for UI and Server
#' @return PCA plot UI
#' @export

pcaModuleUI <- function(id) {
  ns <- NS(id)
  tabPanel("Principal Component Analysis",
           fluidPage(
             fluidRow(
               column(6, plotOutput(ns("pcaPlotOriginal"))),
               column(6, plotOutput(ns("pcaPlotClustering")))
             )
           )
  )
}

#' @title PCA plot Server
#' @description PCA
#' @param id module id for UI and Server
#' @param normCount gene expression matrix
#' @param colData sample information
#' @return PCA plot Server
#' @export

pcaModuleServer <- function(id, normCount, colData) {
  moduleServer(id, function(input, output, session) {
    rownames(normCount) <- normCount$"GeneSymbol"
    normCount <- normCount[, -1]
    #colnames(normCount) <- gsub("\\.", "-", colnames(normCount))
    pcaResult <- prcomp(t(normCount), scale. = TRUE)

    output$pcaPlotOriginal <- renderPlot({
      createPCAPlot(pcaResult, colData, enableClustering = FALSE)
    })
    
    # Render clustering (enableClustering = TRUE)
    output$pcaPlotClustering <- renderPlot({
      createPCAPlot(pcaResult, colData, enableClustering = TRUE)
    })
    
    outputOptions(output, "pcaPlotOriginal", suspendWhenHidden = FALSE)
    outputOptions(output, "pcaPlotClustering", suspendWhenHidden = FALSE)

  })
}