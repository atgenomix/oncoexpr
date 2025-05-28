#' @import shiny
#' @import WGCNA
#' @import DT
#' @import plotly
#' @import sparklyr
NULL

# Utility function: compute Pearson correlation matrix via Spark
sparkCorMatrix <- function(expr_df, sc, transpose = FALSE, feature_col = "features") {
  # expr_df: R data.frame with samples in rows, genes in columns
  #            or if transpose = TRUE, genes in rows, samples in columns
  if (transpose) expr_df <- as.data.frame(t(expr_df))
  tbl <- copy_to(sc, expr_df, name = "expr_tbl", overwrite = TRUE)
  assembler <- ft_vector_assembler(
    tbl,
    input_cols = colnames(expr_df),
    output_col = feature_col
  )
  vec_tbl <- ml_transform(assembler, tbl)
  # Compute Pearson correlation matrix and return as R matrix
  corr_m <- ml_corr(vec_tbl, col = feature_col, method = "pearson")
  corr_mat <- as.matrix(corr_m)  # convert SparkML Matrix to dense R matrix
  # Set appropriate row and column names
  if (transpose) {
    samples <- rownames(expr_df)
    dimnames(corr_mat) <- list(samples, samples)
  } else {
    genes <- colnames(expr_df)
    dimnames(corr_mat) <- list(genes, genes)
  }
  corr_mat
}

#' @title Sample Clustering Shiny Module (Spark-accelerated)
#' @rdname sampleClust
#' @export
sampleClustUI <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("dendro"), height = "400px")
  )
}

#' @rdname sampleClust
#' @export
sampleClustServer <- function(id, exprData, cutHeight) {
  moduleServer(id, function(input, output, session) {
    sampleTree <- reactive({
      df <- exprData()
      sc <- sparklyr::spark_connection()  # establish Spark connection
      # Transpose data to compute sample-to-sample correlation
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
  })
}

#' @title Scale-Free Topology Analysis Shiny Module (Spark compute + local WGCNA)
#' @rdname sft
#' @export
sftUI <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("sftPlot"), height = "400px")
  )
}

#' @rdname sft
#' @export
sftServer <- function(id, exprData, powerRange, rsqCut, selectedPower) {
  moduleServer(id, function(input, output, session) {
    sft_vals <- reactive({
      df <- exprData()
      sc <- sparklyr::spark_connection()  # Spark connection
      # Compute gene correlation matrix via Spark
      corr_mat <- sparkCorMatrix(df, sc, transpose = FALSE)
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
geneModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("blockSelector")),
    plotOutput(ns("geneModulesPlot"), height = "500px")
  )
}

#' @rdname geneModule
#' @export
geneModuleServer <- function(id, exprData, power, deepSplit, minSize, runTrigger) {
  moduleServer(id, function(input, output, session) {
    modulesObj <- eventReactive(runTrigger(), {
      df <- exprData()
      sc <- sparklyr::spark_connection()  # Spark connection
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

    modulesObj
  })
}

#' @title Gene List by Module Color Shiny Module
#' @rdname geneList
#' @export
# (geneListUI/Server unchanged)
