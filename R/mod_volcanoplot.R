#' Interactive Plots UI
#'
#' Creates a user interface component for displaying interactive volcano, scatter, and violin plots.
#'
#' @param id A character string specifying the module ID.
#'
#' @return A Shiny tag list containing the UI layout for interactive plots.
#'
#' @examples
#' \dontrun{
#'   ui <- fluidPage(
#'     interactivePlotsUI("myPlotModule")
#'   )
#' }
#'
#' @export
#' 
#' 

interactivePlotsUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(
      style = "background-color: #f2f2f2; padding: 10px;",
      titlePanel(""),
      fluidRow(
        # Left: Volcano Plot block
        column(width = 8,
               div(style = "background-color: white; border: 2px solid #66CCCC; padding: 10px; margin-bottom: 10px;",
                   h4("Volcano Plot"),
                   
                   girafeOutput(ns("volcanoPlot"), width = "100%", height = "620px")
               )
        ),
        # Right: Scatter Plot (top) and Violin Plot (bottom)
        column(width = 4,
               div(style = "background-color: white; border: 2px solid #66CCCC; padding: 10px; margin-bottom: 10px;",
                   h4("Scatter Plot"),
                   withSpinner(girafeOutput(ns("scatterPlot"), width = "100%", height = "300px"))
               ),
               div(style = "background-color: white; border: 2px solid #66CCCC; padding: 10px;",
                   h4("Violin Plot"),
                   withSpinner(plotOutput(ns("geneViolinPlot"), width = "100%", height = "300px"))
               )
        )
      )
    )
  )
}

#' Interactive Plots Server Module
#'
#' This module creates interactive volcano, scatter, and violin plots based on the provided 
#' differential expression data and expression data in long format. The volcano plot can be adjusted 
#' using external slider parameters passed via \code{params}.
#'
#' @param id A character string specifying the module ID.
#' @param volcanoData A reactive expression that returns a data frame of differential expression results.
#'   The data frame should contain at least the columns \code{GeneSymbol}, \code{PValue}, and \code{logFC}.
#' @param exprData A reactive expression that returns a data frame in long format with columns including 
#'   \code{GeneSymbol}, \code{sample}, \code{expression}, and \code{group}.
#' @param params A reactive expression that returns a list of parameters for the volcano plot. The list should 
#'   include the following elements:
#'   \itemize{
#'     \item \code{lfc_cut} (numeric): Fold Change threshold.
#'     \item \code{pval_cut} (numeric): p-value threshold.
#'     \item \code{pointSize} (numeric): Size of the points.
#'     \item \code{ptAlpha} (numeric): Transparency level of the points.
#'     \item \code{labelSize} (numeric): Font size for gene labels.
#'     \item \code{topN} (numeric): Number of top genes to label.
#'     \item \code{use_adjP} (logical): Whether to use adjusted p-values.
#'   }
#' @param selectedGen selected single gene from DEG filtered table 
#' @return This module does not return a value; instead, it creates server-side outputs for interactive plots.
#'
#' @examples
#' \dontrun{
#'   interactivePlotsServer("myPlotModule",
#'     volcanoData = reactive({ my_deg_data }),
#'     exprData = reactive({ my_long_expr }),
#'     params = reactive({
#'       list(lfc_cut = 1, pval_cut = 0.05, pointSize = 2, ptAlpha = 0.6, labelSize = 3, topN = 20, use_adjP = FALSE)
#'     })
#'   )
#' }
#'
#' @export


interactivePlotsServer <- function(id, volcanoData, exprData, params, selectedGene=NULL) {
  moduleServer(id, function(input, output, session) {
    
    persistent_gene <- reactiveVal(NULL)
    
    observeEvent(input$volcanoPlot_selected, {
      if (!is.null(input$volcanoPlot_selected) && input$volcanoPlot_selected != "")
        persistent_gene(input$volcanoPlot_selected)
    })
    
    observeEvent(input$scatterPlot_selected, {
      if (!is.null(input$scatterPlot_selected) && input$scatterPlot_selected != "")
        persistent_gene(input$scatterPlot_selected)
    })
    
    observeEvent(selectedGene(), {
      if (!is.null(selectedGene()) && length(selectedGene()) > 0) {
        persistent_gene(selectedGene()[1])
      }
    })
    
    output$volcanoPlot <- renderGirafe({
      p <- ggvolcano_custom_interactive(
        df        = volcanoData, 
        geneName  = volcanoData$"GeneSymbol",
        pValCol   = "PValue", 
        logFCCol  = "logFC", 
        lfc_cut   = params()$lfc_cut, 
        pval_cut  = params()$pval_cut, 
        title     = "",
        topN      = params()$topN, 
        pointSize = params()$pointSize, 
        ptAlpha   = params()$ptAlpha, 
        labelSize = params()$labelSize,
        highlight = persistent_gene()
      )
      girafe(ggobj = p,
             options = list(
               opts_zoom(max = 5),
               opts_selection(type = "single", only_shiny = FALSE)
               
             ))
      
    })
    outputOptions(output, "volcanoPlot", suspendWhenHidden = FALSE)
    
    output$scatterPlot <- renderGirafe({
      scatterData_local <- exprData %>%
        group_by(GeneSymbol) %>%
        summarise(
          mean1 = mean(expression[group == "control"]),
          mean2 = mean(expression[group == "case"])
        ) %>% ungroup()
      
      current_gene <- persistent_gene()
      
      if (is.null(current_gene) || current_gene == "") {
        scatterData_local <- scatterData_local %>% mutate(highlight = "normal")
        p <- ggplot(scatterData_local, aes(x = log2(mean1), y = log2(mean2))) +
          geom_point_interactive(aes(tooltip = GeneSymbol, data_id = GeneSymbol),
                                 size = 2, color = "black", alpha = 1) +
          labs(x = "Control (log2 Mean Expression)", y = "Case (log2 Mean Expression)") +
          theme_minimal()
      } else {
        nonhighlight_data <- scatterData_local %>% filter(GeneSymbol != current_gene)
        highlight_data    <- scatterData_local %>% filter(GeneSymbol == current_gene)
        
        p <- ggplot() +
          geom_point_interactive(data = nonhighlight_data,
                                 aes(x = log2(mean1), y = log2(mean2), tooltip = GeneSymbol, data_id = GeneSymbol),
                                 size = 2, color = "grey", alpha = 0.2) +
          geom_point_interactive(data = highlight_data,
                                 aes(x = log2(mean1), y = log2(mean2), tooltip = GeneSymbol, data_id = GeneSymbol),
                                 size = 2, color = "red", alpha = 1) +
          labs(x = "Control (log2 Mean Expression)", y = "Case (log2 Mean Expression)") +
          theme_minimal()
      }
      
      girafe(ggobj = p,
             options = list(
               opts_zoom(max = 5),
               opts_selection(type = "single", only_shiny = FALSE)
             ))
    })
    outputOptions(output, "scatterPlot", suspendWhenHidden = FALSE)

    output$geneViolinPlot <- renderPlot({
      current_gene <- persistent_gene()
      if (is.null(current_gene) || current_gene == "") {
        plot.new()
        title("Please click a gene from the Volcano, Scatter Plot, or Table")
      } else {
        data_selected <- exprData %>% filter(GeneSymbol == current_gene)
        ggplot(data_selected, aes(x = group, y = expression, fill = group)) +
          geom_violin(trim = FALSE, color = NA) +
          geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
          labs(title = paste("Expression of", current_gene),
               x = "Group", y = "Expression Level") +
          theme_minimal() +
          theme(legend.position = "none",
                plot.title = element_text(hjust = 0.5))
      }
    })
    
  })
}
