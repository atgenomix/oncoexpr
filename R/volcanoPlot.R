computeExprData <- function(exprMatrix, groupList) {
  # Check that the length of groupList matches the number of samples
  if(length(groupList) != ncol(exprMatrix)) {
    stop("Length of groupList must equal the number of columns in exprMatrix")
  }
  
  # Convert the matrix to a data frame and add gene names
  df <- as.data.frame(exprMatrix)
  df$gene <- rownames(exprMatrix)
  
  # Convert to long format using tidyr::pivot_longer
  library(tidyr)
  long_df <- pivot_longer(df, cols = -gene, names_to = "sample", values_to = "expression")
  
  # Add group info: assume order of groupList corresponds to the column order of exprMatrix
  long_df$group <- groupList[match(long_df$sample, colnames(exprMatrix))]
  
  return(long_df)
}



interactivePlotsUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Set page background to light gray
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


interactivePlotsServer <- function(id, exactTestData, exprData) {
  moduleServer(id, function(input, output, session) {
    
    # Create a reactiveVal to store the last hovered gene
    persistent_gene <- reactiveVal(NULL)
    
    # Update persistent_gene when hovering over the Volcano Plot
    observeEvent(input$volcanoPlot_hovered, {
      if (!is.null(input$volcanoPlot_hovered) && input$volcanoPlot_hovered != "")
        persistent_gene(input$volcanoPlot_hovered)
    })
    
    # Update persistent_gene when hovering over the Scatter Plot
    observeEvent(input$scatterPlot_hovered, {
      if (!is.null(input$scatterPlot_hovered) && input$scatterPlot_hovered != "")
        persistent_gene(input$scatterPlot_hovered)
    })
    
    # Render Volcano Plot with dynamic color based on persistent_gene
    output$volcanoPlot <- renderGirafe({
      current_gene <- persistent_gene()
      if (is.null(current_gene) || current_gene == "") {
        volcanoData2 <- exactTestData %>% mutate(highlight = "normal")
      } else {
        volcanoData2 <- exactTestData %>% 
          mutate(highlight = ifelse(gene == current_gene, "highlight", "normal"))
      }
      
      p <- ggplot(volcanoData2, aes(x = logFC, y = -log10(pvalue))) +
        geom_point_interactive(aes(tooltip = gene, data_id = gene, color = highlight), size = 4) +
        scale_color_manual(values = c("highlight" = "red", "normal" = "black")) +
        labs(x = "Log Fold Change", y = "-log10(p-value)") +
        theme_minimal()
      
      girafe(ggobj = p,
             options = list(
               opts_zoom(max = 5),
               opts_hover_inv(css = "opacity:0.3;")
             ))
    })
    
    # Render Scatter Plot with dynamic color based on persistent_gene
    output$scatterPlot <- renderGirafe({
      current_gene <- persistent_gene()
      scatterData_local <- exprData %>%
        group_by(gene) %>%
        summarise(
          mean1 = mean(expression[group == "Group1"]),
          mean2 = mean(expression[group == "Group2"])
        ) %>% ungroup()
      
      if (is.null(current_gene) || current_gene == "") {
        scatterData_local <- scatterData_local %>% mutate(highlight = "normal")
      } else {
        scatterData_local <- scatterData_local %>% 
          mutate(highlight = ifelse(gene == current_gene, "highlight", "normal"))
      }
      
      p <- ggplot(scatterData_local, aes(x = mean1, y = mean2)) +
        geom_point_interactive(aes(tooltip = gene, data_id = gene, color = highlight), size = 4) +
        scale_color_manual(values = c("highlight" = "red", "normal" = "black")) +
        labs(x = "Group1 Mean Expression", y = "Group2 Mean Expression") +
        theme_minimal()
      
      girafe(ggobj = p,
             options = list(
               opts_zoom(max = 5),
               opts_hover_inv(css = "opacity:0.3;")
             ))
    })
    
    # Render Violin Plot based on persistent_gene
    output$geneViolinPlot <- renderPlot({
      selected_gene <- persistent_gene()
      if (is.null(selected_gene) || selected_gene == "") {
        plot.new()
        title("Please hover over a gene from the Volcano or Scatter Plot")
      } else {
        data_selected <- exprData %>% filter(gene == selected_gene)
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





# Example expression matrix
set.seed(123)
exprMatrix <- matrix(rnorm(50), nrow = 10, ncol = 5)
rownames(exprMatrix) <- paste0("Gene", 1:10)
colnames(exprMatrix) <- paste0("Sample", 1:5)

# Example group assignments for each sample
groupList <- c("Group1", "Group1", "Group2", "Group2", "Group1")

# Convert to long format data frame
exprData <- computeExprData(exprMatrix, groupList)
print(exprData)


# set.seed(123)

# # Generate Volcano Plot data: genes with log fold change and p-value
# volcanoData <- data.frame(
#   gene = paste0("Gene", 1:10),
#   logFC = rnorm(10, 0, 2),
#   pvalue = runif(10, 0.001, 0.1)
# )

# # Generate expression data for each gene (3 samples per group)
# exprData <- do.call(rbind, lapply(1:nrow(volcanoData), function(i) {
#   gene_name <- volcanoData$gene[i]
#   effect <- volcanoData$logFC[i]
#   group1_mean <- 10 - effect/2
#   group2_mean <- 10 + effect/2
#   data.frame(
#     gene = gene_name,
#     group = rep(c("Group1", "Group2"), each = 3),
#     expression = c(rnorm(3, mean = group1_mean, sd = 2),
#                    rnorm(3, mean = group2_mean, sd = 2))
#   )
# }))

ui <- fluidPage(
  interactivePlotsUI("plotModule")
)

server <- function(input, output, session) {
  interactivePlotsServer("plotModule", volcanoData = volcanoData, exprData = exprData)
}

shinyApp(ui = ui, server = server)
