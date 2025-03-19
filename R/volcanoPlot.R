library(shiny)
library(ggplot2)
library(ggiraph)
library(dplyr)
library(tidyr)

computeExprData <- function(exprMatrix = normCount, colData = colData) {
  df <- as.data.frame(exprMatrix)
  colData$mainCode <- colnames(df)[-which(colnames(df)=="GeneSymbol")]
  print(colnames(df))
  long_df <- pivot_longer(df, cols = -GeneSymbol , names_to = "sample", values_to = "expression")
  long_df$group <- colData$subCode[match(long_df$sample, colData$mainCode)]
  
  return(long_df)
}


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


interactivePlotsServer <- function(id, volcanoData, exprData, params) {
  moduleServer(id, function(input, output, session) {
    
    # 其他 reactive 與 persistent_gene 的設定…
    persistent_gene <- reactiveVal(NULL)
    observeEvent(input$volcanoPlot_hovered, {
      if (!is.null(input$volcanoPlot_hovered) && input$volcanoPlot_hovered != "")
        persistent_gene(input$volcanoPlot_hovered)
    })
    observeEvent(input$scatterPlot_hovered, {
      if (!is.null(input$scatterPlot_hovered) && input$scatterPlot_hovered != "")
        persistent_gene(input$scatterPlot_hovered)
    })
    
    # 使用傳入的 params 更新 Volcano Plot
    output$volcanoPlot <- renderGirafe({
      # 讀取外部 slider 參數 (注意 params 為 reactive)
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
               opts_hover_inv(css = "opacity:0.3;")
             ))
    })
    
    # 以下 scatterPlot 與 geneViolinPlot 部分不受 slider 影響，可照原來邏輯處理
    output$scatterPlot <- renderGirafe({
        current_gene <- persistent_gene()
        scatterData_local <- exprData %>%
        group_by(GeneSymbol) %>%
        summarise(
            mean1 = mean(expression[group == "control"]),
            mean2 = mean(expression[group == "case"])
        ) %>% ungroup()
      
        if (is.null(current_gene) || current_gene == "") {
            scatterData_local <- scatterData_local %>% mutate(highlight = "normal")
        } else {
            scatterData_local <- scatterData_local %>% 
            mutate(highlight = ifelse(GeneSymbol == current_gene, "highlight", "other"))

        }

        scatterData_local <- scatterData_local %>% 
            mutate(highlight = factor(highlight, levels = c("highlight", "other", "normal")))
        p <- ggplot(scatterData_local, aes(x = mean1, y = mean2)) +
            geom_point_interactive(aes(tooltip = GeneSymbol, data_id = GeneSymbol, color = highlight, alpha = factor(highlight)), size = 2) +
            scale_color_manual(values = c("highlight" = "red", "other" = "grey", "normal" = "black"), drop = FALSE) +
            scale_alpha_manual(values = c("highlight" = 1, "other" = 0.2, "normal" = 1)) +
            labs(x = "Group1 Mean Expression", y = "Group2 Mean Expression") +
            theme_minimal()
        
        girafe(ggobj = p,
                options = list(
                opts_zoom(max = 5),
                opts_hover_inv(css = "opacity:0.3;")
                ))
    })
    
    output$geneViolinPlot <- renderPlot({
      selected_gene <- persistent_gene()
      if (is.null(selected_gene) || selected_gene == "") {
        plot.new()
        title("Please hover over a gene from the Volcano or Scatter Plot")
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


normCount <- read.csv("/Users/charleschuang/Desktop/NormCountSalmon.csv")
exactTest <- read.csv("/Users/charleschuang/Desktop/ExactTestSalmon.csv",row.names=1)
colnames(normCount)[colnames(normCount) == "genes"] <- "GeneSymbol"
colnames(exactTest)[colnames(exactTest) == "genes"] <- "GeneSymbol"

volcanoData <- exactTest

classification <- c(rep("control",3), rep("case",3))
colData <- data.frame("mainCode" = colnames(normCount)[-1], "subCode" = classification)
rownames(colData) <- colnames(normCount)[-1]



#expr_mat <- as.matrix(normCount[, -which(colnames(normCount) == "GeneSymbol")])
exprData <- computeExprData( normCount, colData)
View(exprData)
View(volcanoData)
# ---- Main App ----
ui <- fluidPage(
  tabPanel(
        title = "Differential Expression Analysis",
        sidebarLayout(
          sidebarPanel(
            sliderInput("lfc_cut", "Fold Change Threshold (log2):", 
                        min = 0, max = 10, value = 1, step = 0.1),
            sliderInput("pval_cut", "p-value Threshold:", 
                        min = 0.001, max = 1, value = 0.05, step = 0.001),
            sliderInput("pointSize", "Point Size:", 
                        min = 1, max = 5, value = 2, step = 0.5),
            sliderInput("ptAlpha", "Transparent:", 
                        min = 0.1, max = 1, value = 0.6, step = 0.1),
            sliderInput("labelSize", "Gene Label Size:", 
                        min = 1, max = 6, value = 3, step = 0.5),
            numericInput("topN", "Label N number of Genes (0 is no labeling):", 
                         value = 100, min = 0, max = 1000),
            checkboxInput("use_adjP", "Use Adjusted P-value?", value = FALSE),
            actionButton("run_DEG", "Run DEG analysis")
          ),
          mainPanel(
            tabsetPanel(
              tabPanel("Volcano Plot interaction",
                      #  plotOutput("volcano_plot", height = "600px")
                      interactivePlotsUI("plotVolcano")
              ),
              tabPanel("DEG Table", 
                       DT::dataTableOutput('DEG_table', width = "100%")
              )
            )
          )
        )
      )
)

server <- function(input, output, session) {
  # 假設 volcanoData 與 exprData 已正確定義
    params <- reactive({
      list(
        lfc_cut    = input$lfc_cut,
        pval_cut   = input$pval_cut,
        pointSize  = input$pointSize,
        ptAlpha    = input$ptAlpha,
        labelSize  = input$labelSize,
        topN       = input$topN,
        use_adjP   = input$use_adjP
      )
    })
  interactivePlotsServer("plotVolcano", 
                           volcanoData = volcanoData, 
                           exprData = exprData, params)
}
shinyApp(ui = ui, server = server)



# Example expression matrix
set.seed(123)
exprMatrix <- matrix(rnorm(50), nrow = 10, ncol = 5)
rownames(exprMatrix) <- paste0("Gene", 1:10)
colnames(exprMatrix) <- paste0("Sample", 1:5)

# Example group assignments for each sample
groupList <- c("Group1", "Group1", "Group2", "Group2", "Group1")

# Convert to long format data frame
exprData <- computeExprData(exprMatrix, colData)
print(exprData)


# set.seed(123)

# Generate Volcano Plot data: genes with log fold change and p-value
volcanoData <- data.frame(
  gene = paste0("Gene", 1:10),
  logFC = rnorm(10, 0, 2),
  pvalue = runif(10, 0.001, 0.1)
)

# Generate expression data for each gene (3 samples per group)
exprData <- do.call(rbind, lapply(1:nrow(volcanoData), function(i) {
  gene_name <- volcanoData$gene[i]
  effect <- volcanoData$logFC[i]
  group1_mean <- 10 - effect/2
  group2_mean <- 10 + effect/2
  data.frame(
    gene = gene_name,
    group = rep(c("Group1", "Group2"), each = 3),
    expression = c(rnorm(3, mean = group1_mean, sd = 2),
                   rnorm(3, mean = group2_mean, sd = 2))
  )
}))


head(volcanoData)

p <- ggvolcano_custom(
        df        = volcanoData, 
        geneName  = volcanoData$"GeneSymbol",
        pValCol   = "PValue", 
        logFCCol  = "logFC", 
        lfc_cut   = 1, 
        pval_cut  = 0.05, 
        title     = "Volcano Plot",
        topN      = 100, 
        pointSize = 3, 
        ptAlpha   = 0.5, 
        labelSize = 5
      )
girafe(ggobj = p,
        options = list(
        opts_zoom(max = 5),
        opts_hover_inv(css = "opacity:0.3;")
        ))