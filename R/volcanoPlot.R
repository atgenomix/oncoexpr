#' Compute Expression Data in Long Format
#'
#' This function converts a given expression matrix and associated sample metadata 
#' (colData) into a long-format data frame suitable for downstream plotting.
#' The function extracts gene expression values from the matrix and assigns each 
#' sample a group based on the matching of sample names with the metadata.
#'
#' @param exprMatrix A numeric matrix or data frame of expression data, where rows represent genes
#'   and columns represent samples. It is expected that one column is named "GeneSymbol" (or that row names
#'   contain gene symbols) which will be used to annotate the genes.
#' @param colData A data frame containing sample metadata. It must include a column named 
#'   \code{subCode} representing the group labels. This function will also add a \code{mainCode} column 
#'   to \code{colData} based on the column names of \code{exprMatrix} (excluding the "GeneSymbol" column).
#'
#' @return A data frame in long format with the following columns:
#' \describe{
#'   \item{GeneSymbol}{Gene names.}
#'   \item{sample}{Sample identifiers.}
#'   \item{expression}{Gene expression values.}
#'   \item{group}{Group label corresponding to each sample, taken from \code{colData$subCode}.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Assuming normCount is an expression matrix and colData is a data frame with a 'subCode' column:
#'   long_expr <- transfExprFormat(normCount, colData)
#' }
#'
#' @export
transfExprFormat <- function(exprMatrix = normCount, colData = colData) {
  df <- as.data.frame(exprMatrix)
  colData$mainCode <- colnames(df)[-which(colnames(df) == "GeneSymbol")]
  print(colnames(df))
  long_df <- tidyr::pivot_longer(df, cols = -GeneSymbol, names_to = "sample", values_to = "expression")
  long_df$group <- colData$subCode[match(long_df$sample, colData$mainCode)]
  
  return(long_df)
}

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
      titlePanel("Interactive Volcano, Scatter, and Violin Plots"),
      fluidRow(
        # Left: Volcano Plot block
        column(width = 8,
               div(style = "background-color: white; border: 2px solid #66CCCC; padding: 10px; margin-bottom: 10px;",
                   h4("Volcano Plot"),
                   
                   girafeOutput(ns("volcanoPlot"), width = "100%", height = "600px")
               )
        ),
        # Right: Scatter Plot (top) and Violin Plot (bottom)
        column(width = 4,
               div(style = "background-color: white; border: 2px solid #66CCCC; padding: 10px; margin-bottom: 10px;",
                   h4("Scatter Plot"),
                   girafeOutput(ns("scatterPlot"), width = "100%", height = "300px")
               ),
               div(style = "background-color: white; border: 2px solid #66CCCC; padding: 10px;",
                   h4("Violin Plot"),
                   plotOutput(ns("geneViolinPlot"), width = "100%", height = "300px")
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
#'
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


interactivePlotsServer <- function(id, volcanoData, exprData, params) {
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
    

    output$volcanoPlot <- renderGirafe({
      p <- ggvolcano_custom_interactive(
        df        = volcanoData, 
        geneName  = volcanoData$"GeneSymbol",
        pValCol   = "PValue", 
        logFCCol  = "logFC", 
        lfc_cut   = params()$lfc_cut, 
        pval_cut  = params()$pval_cut, 
        title     = "Volcano Plot",
        topN      = params()$topN, 
        pointSize = params()$pointSize, 
        ptAlpha   = params()$ptAlpha, 
        labelSize = params()$labelSize
      )
      girafe(ggobj = p,
             options = list(
               opts_zoom(max = 5),
               opts_selection(type = "single", only_shiny = FALSE)
             ))
    })
    
   output$scatterPlot <- renderGirafe({
        current_gene <- persistent_gene()
        scatterData_local <- exprData %>%
            group_by(GeneSymbol) %>%
            summarise(
            mean1 = mean(expression[group == "control"]),
            mean2 = mean(expression[group == "case"])
            ) %>% ungroup()
        
        # 如果沒有選取基因，直接全部以 normal 繪製
        if (is.null(current_gene) || current_gene == "") {
            scatterData_local <- scatterData_local %>% mutate(highlight = "normal")
            p <- ggplot(scatterData_local, aes(x = log2(mean1), y = log2(mean2))) +
            geom_point_interactive(aes(tooltip = GeneSymbol, data_id = GeneSymbol),
                                    size = 2, color = "black", alpha = 1) +
            labs(x = "Group1 Mean Expression", y = "Group2 Mean Expression") +
            theme_minimal()
        } else {
            # 有選取基因時，分成兩層：非選取與選取
            nonhighlight_data <- scatterData_local %>% filter(GeneSymbol != current_gene)
            highlight_data    <- scatterData_local %>% filter(GeneSymbol == current_gene)
            
            p <- ggplot() +
            geom_point_interactive(data = nonhighlight_data,
                                    aes(x = mean1, y = mean2, tooltip = GeneSymbol, data_id = GeneSymbol),
                                    size = 2, color = "grey", alpha = 0.2) +
            geom_point_interactive(data = highlight_data,
                                    aes(x = mean1, y = mean2, tooltip = GeneSymbol, data_id = GeneSymbol),
                                    size = 2, color = "red", alpha = 1) +
            labs(x = "Control Mean Expression", y = "Case Mean Expression") +
            theme_minimal()
        }
        
        girafe(ggobj = p,
                options = list(
                opts_zoom(max = 5),
                opts_selection(type = "single", only_shiny = FALSE)
                ))
    })

    
    output$geneViolinPlot <- renderPlot({
      selected_gene <- persistent_gene()
      if (is.null(selected_gene) || selected_gene == "") {
        plot.new()
        title("Please click a gene from the Volcano or Scatter Plot")
      } else {
        data_selected <- exprData %>% filter(GeneSymbol == selected_gene)
        ggplot(data_selected, aes(x = group, y = expression, fill = group)) +
          geom_violin(trim = FALSE, color = NA) +
          geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
          labs(title = paste("Expression of", selected_gene),
               x = "Group", y = "Expression Level") +
          theme_minimal() +
          theme(legend.position = "none",
                plot.title = element_text(hjust = 0.5))
      }
    })
    
  })
}
