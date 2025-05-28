#' @import shiny
#' @import WGCNA
#' @import DT
#' @import plotly
NULL

#' @title Sample Clustering Shiny Module
#' @description
#' The `sampleClustUI()` function generates the UI for displaying a sample dendrogram,
#' and `sampleClustServer()` computes pairwise distances, performs hierarchical clustering,
#' and renders a dendrogram with a user-defined cutoff line.
#'
#' @param id `character(1)` Module namespace identifier.
#' @param exprData `reactive` returning a data frame of numeric expression data
#'   (samples in rows, genes in columns).
#' @param distMethod `reactive` specifying the distance metric for `dist()`.
#' @param cutHeight `reactive` numeric value indicating the height at which to draw
#'   the horizontal cutoff line on the dendrogram.
#'
#' @return
#' - `sampleClustUI()`: A Shiny UI definition (`tagList`).
#' - `sampleClustServer()`: A server-side function rendering the dendrogram.
#'
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
sampleClustServer <- function(id, exprData, distMethod, cutHeight) {
  moduleServer(id, function(input, output, session) {
    sampleTree <- reactive({
      dist(exprData(), method = distMethod()) |> hclust(method = "average")
    })
    maxHeight <- reactive({
      max(sampleTree()$height)
    })
    sampleClusters <- reactive({
      tree <- sampleTree()
      cutreeStatic(tree, cutHeight = cutHeight(), minSize = 1)
    })
    
    filteredExprData <- reactive({
      clusters <- sampleClusters()
      message("Current cutHeight: ", cutHeight())
      message("sample clusters: ", clusters)      
      keepGroup <- which.max(table(clusters))
      exprData()[clusters == keepGroup, , drop = FALSE]
      
    })
    message("filteredExprData: ", dim(filteredExprData))      
    output$dendro <- renderPlot({
      tree <- sampleTree()
      plot(tree,
           main = "Sample Clustering to Detect Outliers",
           sub = "", xlab = "",
           cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.5)
      abline(h = cutHeight(), col = "red", lwd = 2)
    })

    return(list(
      filteredExpr = filteredExprData,
      maxHeight = maxHeight
    ))

  })
}

#' @title Scale-Free Topology Analysis Shiny Module
#' @description
#' The `sftUI()` function generates UI for plotting scale-free topology fit indices,
#' and `sftServer()` performs `WGCNA::pickSoftThreshold()` over a range of powers,
#' visualizing model fit and highlighting the selected power and R^2 cutoff.
#'
#' @param id `character(1)` Module namespace identifier.
#' @param exprData `reactive` returning expression data as a numeric data frame.
#' @param powerRange `reactive` numeric vector of length 2 specifying the min and max
#'   soft-thresholding powers to evaluate.
#' @param rsqCut `reactive` numeric value for the R^2 cutoff used in `pickSoftThreshold()`.
#' @param selectedPower `reactive` numeric indicating the chosen power to highlight on the plot.
#'
#' @return
#' - `sftUI()`: A Shiny UI definition (`tagList`).
#' - `sftServer()`: A server-side function rendering the scale-free fit plot.
#'
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
    sft <- reactive({
      withProgress(message = "Calculating scale-free fit...", value = 0, {
        pickSoftThreshold(
          exprData(),
          powerVector = seq(powerRange()[1], powerRange()[2]),
          RsquaredCut = rsqCut(),
          verbose = 0
        )
      })
    })
    output$sftPlot <- renderPlot({
      fit <- sft()$fitIndices
      plot(fit[, 1], fit[, 2],
           xlab = "Soft Threshold (power)",
           ylab = "Scale Free Topology Model Fit, signed R^2",
           type = "n",
           main = "Scale Independence")
      text(fit[, 1], fit[, 2], labels = fit[, 1], col = "red")
      abline(v = selectedPower(), col = "blue", lwd = 2)
      abline(h = rsqCut(), col = "red", lwd = 2)
    })
  })
}

#' @title Module Detection and Dendrogram Shiny Module
#' @description
#' `geneModuleUI()` displays a gene dendrogram and module color strips,
#' and `geneModuleServer()` calls `WGCNA::blockwiseModules()` to detect co-expression
#' modules, rendering the dendrogram with colored modules. Returns the modules object
#' for downstream use.
#'
#' @param id `character(1)` Module namespace identifier.
#' @param exprData `reactive` numeric expression data.
#' @param power, deepSplit, minSize `reactive` parameters passed to `blockwiseModules()`.
#' @param runTrigger `reactive` trigger (e.g., actionButton) to re-run module detection.
#'
#' @return
#' - `geneModuleServer()`: A `reactive` blockwiseModules result object.
#'
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
      #nGenes <- ncol(exprData())
      withProgress(message = "Detecting modules...", value = 0, {
        blockwiseModules(
          exprData(),
          power = power(),
          corFnc = WGCNA::cor,
          deepSplit = deepSplit(),
          minModuleSize = minSize(),
          numericLabels = FALSE,
          maxBlockSize = 1000,
          verbose = 0
        )
      })
    })
    ########
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
    ########
    output$geneModulesPlot <- renderPlot({
      obj <- modulesObj()
      validate(need(!is.null(obj), "Click 'Run WGCNA' to detect modules."))
      req(obj, input$block)
      b_ <- as.integer(input$block)
      print(paste("Block:", b_))
      tree <- obj$dendrograms[[b_]] # Use the first dendrogram
      colors <- labels2colors(obj$colors)
      names(colors) <- colnames(exprData())
      subcolors <- obj$colors[obj$blockGenes[[b_]]]
      message("Number of modules detected: ", length(unique(colors)))
      message("Total genes colored: ", length(colors))
      message("exprData dimensions: ", 
              paste(nrow(exprData()), ncol(exprData()), sep = " x "))
      message("numbers of blocks: ",  length(obj$dendrograms))
      message("number of genes in the first block: ", tree$size)
      message("number of subcolors: ", length(subcolors) )
      message("number of genes in the first block: ", length(obj$blockGenes[[b_]]))
      plotDendroAndColors(
        tree, subcolors, "Module",
        dendroLabels = FALSE,
        hang = 0.03,
        addGuide = TRUE,
        guideHang = 0.05,
        main = paste("Gene Dendrogram — Block", b_)
      )
    })
    modulesObj
  })
}

#' @title Gene List by Module Color Shiny Module
#' @description
#' `geneListUI()` provides a UI for selecting a module color and displaying
#' the corresponding gene list; `geneListServer()` populates the table and
#' allows downloading the gene list as CSV.
#'
#' @param id `character(1)` Module namespace identifier.
#' @param exprData `reactive` numeric expression data.
#' @param modulesObj `reactive` blockwiseModules result from `geneModuleServer()`.
#'
#' @return
#' - UI: A color selector, a data table, and a download button.
#' - Server: Renders the DataTable and handles file download.
#'
#' @rdname geneList
#' @export
geneListUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("moduleSelector")),
    DTOutput(ns("geneTable")),
    fluidRow(
      column(6, downloadButton(ns("downloadGenes"), "Download Selected Module")),
      column(6, downloadButton(ns("downloadAll"),   "Download All Modules"))
    )
    
  )
}

#' @rdname geneList
#' @export
geneListServer <- function(id, exprData, modulesObj) {
  moduleServer(id, function(input, output, session) {
    df_all <- reactive({
      modobj <- modulesObj()
      validate(need(!is.null(modobj), "Detect modules first."))
      genes <- colnames(exprData())
      colors <- labels2colors(modobj$colors)
      data.frame(Gene = genes, ModuleColor = colors, stringsAsFactors = FALSE)
    })
    output$moduleSelector <- renderUI({
      df <- df_all()
      cols <- sort(unique(df$ModuleColor))
      selectInput(session$ns("moduleColor"), "Select Module Color:", choices = cols)
    })
    output$geneTable <- renderDT(
      {
        req(input$moduleColor)
        subset(df_all(), ModuleColor == input$moduleColor)
      },
      options = list(pageLength = 10, lengthChange = FALSE)
    )
    output$downloadGenes <- downloadHandler(
      filename = function() paste0("genes_module_", input$moduleColor, ".csv"),
      content = function(file) {
        write.csv(subset(df_all(), ModuleColor == input$moduleColor), file, row.names = FALSE)
      }
    )
    output$downloadAll <- downloadHandler(
      filename = function() {
        "all_modules.csv"
      },
      content = function(file) {
        write.csv(df_all(), file, row.names = FALSE)
      }
    )
  })
}
