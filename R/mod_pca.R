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
    normCount <- normCount[, setdiff(names(normCount), "GeneSymbol")]
    sample_n <- ncol(normCount)
    if (sample_n < 3) {
      output$pcaPlotOriginal <- renderPlot({
        plot.new()
        text(0.5, 0.5, "Insufficient samples to perform PCA.")
        box(col = "cyan", lwd = 1)
      })
      output$pcaPlotClustering <- renderPlot({
        plot.new()
        text(0.5, 0.5, "Insufficient samples to perform PCA.")
        box(col = "cyan", lwd = 1)
      })
      outputOptions(output, "pcaPlotOriginal",     suspendWhenHidden = FALSE)
      outputOptions(output, "pcaPlotClustering",   suspendWhenHidden = FALSE)
      return()
    }
    #colnames(normCount) <- gsub("\\.", "-", colnames(normCount))
    pcaResult <- prcomp(t(normCount), scale. = TRUE)

    output$pcaPlotOriginal <- renderPlot({
      createPCAPlot(pcaResult, colData, enableClustering = FALSE)+
        theme(
          panel.border = element_rect(color = "cyan", fill = NA, size = 1)
        )
    })
    
    # Render clustering (enableClustering = TRUE)
    output$pcaPlotClustering <- renderPlot({
      createPCAPlot(pcaResult, colData, enableClustering = TRUE)+
        theme(
          panel.border = element_rect(color = "cyan", fill = NA, size = 1)
        )
    })
    
    outputOptions(output, "pcaPlotOriginal", suspendWhenHidden = FALSE)
    outputOptions(output, "pcaPlotClustering", suspendWhenHidden = FALSE)

  })
}