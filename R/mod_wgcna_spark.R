#' @import shiny
#' @import WGCNA
#' @import DT
#' @import plotly
#' @import sparklyr
NULL

# Utility function: compute Pearson correlation matrix via Spark
sparkCorMatrix <- function(expr_df, sc, transpose = FALSE, feature_col = "features") {
  # If a reactive or function is mistakenly passed, evaluate it
  if (is.function(expr_df)) {
    expr_df <- expr_df()
  }
  # Ensure expr_df is a matrix, then convert to data frame after optional transpose
  mat <- as.matrix(expr_df)
  if (transpose) {
    mat <- t(mat)
  }
  expr_df2 <- as.data.frame(mat)

  # Copy to Spark and assemble vector
  tbl <- copy_to(sc, expr_df2, name = "expr_tbl", overwrite = TRUE)
  vec_tbl <- tbl %>%
    ft_vector_assembler(input_cols = colnames(expr_df2), output_col = feature_col)

  # Compute correlation
  corr_m  <- ml_corr(vec_tbl, col = feature_col, method = "pearson")
  corr_mat <- as.matrix(corr_m)

  # Assign appropriate dimnames based on original expr_df
  if (transpose) {
    sample_names <- rownames(expr_df)
    dimnames(corr_mat) <- list(sample_names, sample_names)
  } else {
    gene_names <- colnames(expr_df)
    dimnames(corr_mat) <- list(gene_names, gene_names)
  }
  corr_mat
}

#' @title Sample Clustering Shiny Module (Spark-accelerated)
#' @rdname sampleClust
#' @export
sparksampleClustUI <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("dendro"), height = "400px")
  )
}

#' @rdname sampleClust
#' @export
sparksampleClustServer <- function(id, exprData, cutHeight) {
  moduleServer(id, function(input, output, session) {
    sampleTree <- reactive({
      df <- exprData()
      #sc <- sparklyr::spark_connection()  # establish Spark connection
      # Transpose data to compute sample-to-sample correlation
      message("Computing sample correlation matrix via Spark...")
      corr_mat <- sparkCorMatrix(df, sc, transpose = TRUE)
      dist_mat <- as.dist(1 - corr_mat)
      hclust(dist_mat, method = "average")
    })

    # Render dendrogram plot
    output$dendro <- renderPlot({
      tree <- sampleTree()
      plot(
        tree,
        main = "Sample Clustering (Spark-accelerated)",
        sub = "", xlab = "",
        cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.5
      )
      abline(h = cutHeight(), col = "red", lwd = 2)
    })

    maxHeight <- reactive({
      tree <- sampleTree()
      max(tree$height)
    })

    filteredExpr <- reactive({
      exprData()
    })

    return(list(
      maxHeight    = maxHeight,
      filteredExpr = filteredExpr
    ))

  })
}

#' @title Scale-Free Topology Analysis Shiny Module (Spark compute + local WGCNA)
#' @rdname sft
#' @export
sparksftUI <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("sftPlot"), height = "400px")
  )
}

#' @rdname sft
#' @export
sparksftServer <- function(id, exprData, powerRange, rsqCut, selectedPower) {
  moduleServer(id, function(input, output, session) {
    sft_vals <- reactive({
      df <- exprData()
      #sc <- sparklyr::spark_connection()  # Spark connection
      # Compute gene correlation matrix via Spark
      message("Computing correlation matrix via Spark...")
      corr_mat <- sparkCorMatrix(df, sc, transpose = FALSE)
      message("Correlation matrix computed. Evaluating soft-threshold fit...")
      # Evaluate soft-threshold fit using WGCNA locally
      pickSoftThreshold(
        corr_mat,
        powerVector = seq(powerRange()[1], powerRange()[2]),
        RsquaredCut = rsqCut(),
        verbose = 0,
        dataIsExpr = FALSE  # input is correlation matrix, not raw expression
      )
    })

    output$sftPlot <- renderPlot({
      fit <- sft_vals()$fitIndices
      plot(
        fit[, 1], fit[, 2],
        xlab = "Soft Threshold (power)",
        ylab = "Scale Free Topology Model Fit, signed R^2",
        type = "n",
        main = "Scale Independence"
      )
      text(fit[, 1], fit[, 2], labels = fit[, 1], col = "red")
      abline(v = selectedPower(), col = "blue", lwd = 2)
      abline(h = rsqCut(), col = "red", lwd = 2)
    })
  })
}

#' @title Module Detection Shiny Module (Spark compute correlation + local WGCNA)
#' @rdname geneModule
#' @export
sparkgeneModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("blockSelector")),
    plotOutput(ns("geneModulesPlot"), height = "500px")
  )
}

#' @rdname geneModule
#' @export
sparkgeneModuleServer <- function(id, exprData, power, deepSplit, minSize, runTrigger) {
  moduleServer(id, function(input, output, session) {
    modulesObj <- eventReactive(runTrigger(), {
      df <- exprData()
      #sc <- sparklyr::spark_connection()  # Spark connection
      # Compute gene correlation matrix via Spark
      corr_mat <- sparkCorMatrix(df, sc, transpose = FALSE)
      # Detect modules using blockwiseModules with correlation matrix input
      blockwiseModules(
        corr_mat,
        power = power(),
        corOptions = list(dataIsExpr = FALSE),
        corFnc = WGCNA::cor,
        deepSplit = deepSplit(),
        minModuleSize = minSize(),
        numericLabels = FALSE,
        maxBlockSize = 5000,
        verbose = 0
      )
    })

    output$blockSelector <- renderUI({
      mod <- modulesObj()
      req(mod)
      nblocks <- length(mod$dendrograms)
      selectInput(
        session$ns("block"),
        label = "Select block:",
        choices = setNames(seq_len(nblocks), paste0("Block ", seq_len(nblocks))),
        selected = 1
      )
    })

    output$geneModulesPlot <- renderPlot({
      obj <- modulesObj()
      req(obj, input$block)
      block_index <- as.integer(input$block)
      tree <- obj$dendrograms[[block_index]]
      colors <- labels2colors(obj$colors)
      subcols <- obj$colors[obj$blockGenes[[block_index]]]
      plotDendroAndColors(
        tree, subcols, "Module",
        dendroLabels = FALSE,
        hang = 0.03,
        addGuide = TRUE,
        guideHang = 0.05,
        main = paste("Gene Modules — Block", block_index)
      )
    })

    return(modulesObj)
  })
}

#' @title Gene List by Module Color Shiny Module
#' @rdname geneList
#' @export
# (geneListUI/Server unchanged)




"====="

library(shiny)

#test 
sc <- spark_connect(master = "local")

# When the Shiny session stops, close the connection gracefully
onStop(function() {
  try(spark_disconnect(sc), silent = TRUE)
})

ui <- fluidPage(
  titlePanel("Spark WGCNA ShinyApp"),
  sidebarLayout(
    sidebarPanel(
      numericInput("cutHeight", "Cut Height (Dendrogram):", 0.5, min = 0, step = 0.1),
      numericInput("power", "Soft Threshold Power:", 6, min = 1, max = 20),
      numericInput("deepSplit", "Deep Split:", 2, min = 0, max = 4),
      numericInput("minSize", "Min Module Size:", 30, min = 10),
      actionButton("run", "Run Module Detection")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Sample Clustering", sparksampleClustUI("sampleClust")),
        tabPanel("Scale-Free Analysis", sparksftUI("sft")),
        tabPanel("Gene Module Detection", sparkgeneModuleUI("geneMod"))
      )
    )
  )
)

# Shiny Server
server <- function(input, output, session) {
  # Example expression data
  exprData <- reactive({
    matrix(rnorm(1000), nrow = 20, dimnames = list(paste0("Sample", 1:20), paste0("Gene", 1:50)))
  })

  # Sample clustering module
  sparksampleClustServer("sampleClust", exprData, reactive(input$cutHeight))

  # Scale-free topology module
  sparksftServer("sft", exprData, reactive(c(1, 10)), reactive(0.8), reactive(input$power))

  # Gene module detection module
  sparkgeneModuleServer("geneMod", exprData,
                        reactive(input$power), reactive(input$deepSplit), reactive(input$minSize),
                        reactive(input$run))
}

shinyApp(ui, server)




