# app.R: WGCNA Shiny App with Shiny Modules (Corrected geneListServer)
# ================================
library(shiny)
library(WGCNA)
library(DT)

# Enable WGCNA multithreading and global options
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

#----- Module 1: Sample Clustering --------------------------------
sampleClustUI <- function(id) {
    ns <- NS(id)
    tagList(
        plotOutput(ns("dendro"), height = "400px")
    )
}

sampleClustServer <- function(id, exprData, distMethod, cutHeight) {
    moduleServer(id, function(input, output, session) {
        sampleTree <- reactive({
            dist(exprData(), method = distMethod()) |> hclust(method = "average")
        })
        output$dendro <- renderPlot({
            tree <- sampleTree()
            plot(tree,
                main = "Sample clustering to detect outliers",
                sub = "", xlab = "",
                cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.5
            )
            abline(h = cutHeight(), col = "red", lwd = 2)
        })
    })
}

#----- Module 2: Scale-Free Topology -------------------------------
sftUI <- function(id) {
    ns <- NS(id)
    tagList(
        plotOutput(ns("sftPlot"), height = "400px")
    )
}

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
                main = "Scale independence"
            )
            text(fit[, 1], fit[, 2], labels = fit[, 1], col = "red")
            abline(v = selectedPower(), col = "blue", lwd = 2)
            abline(h = rsqCut(), col = "red", lwd = 2)
        })
    })
}

#----- Module 3: Module Detection & Dendrogram ---------------------
geneModuleUI <- function(id) {
    ns <- NS(id)
    tagList(
        plotOutput(ns("geneModulesPlot"), height = "500px")
    )
}

geneModuleServer <- function(id, exprData, power, deepSplit, minSize, runTrigger) {
    moduleServer(id, function(input, output, session) {
        modulesObj <- eventReactive(runTrigger(), {
            withProgress(message = "Detecting modules...", value = 0, {
                blockwiseModules(
                    exprData(),
                    power = power(),
                    corFnc = WGCNA::cor,
                    deepSplit = deepSplit(),
                    minModuleSize = minSize(),
                    numericLabels = FALSE,
                    verbose = 0
                )
            })
        })
        output$geneModulesPlot <- renderPlot({
            obj <- modulesObj()
            validate(need(!is.null(obj), "Click 'Detect Modules' first"))
            tree <- obj$dendrograms[[1]]
            colors <- labels2colors(obj$colors)
            plotDendroAndColors(
                tree, colors, "Module",
                dendroLabels = FALSE,
                hang = 0.03,
                addGuide = TRUE,
                guideHang = 0.05,
                main = "Gene dendrogram and module colors"
            )
        })
        return(modulesObj)
    })
}

#----- Module 4: Gene List per Module -------------------------------
geneListUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("moduleSelector")),
        DTOutput(ns("geneTable")),
        downloadButton(ns("downloadGenes"), "Download CSV")
    )
}

geneListServer <- function(id, exprData, modulesObj) {
    moduleServer(id, function(input, output, session) {
        df_all <- reactive({
            modobj <- modulesObj()
            validate(need(!is.null(modobj), "Detect modules first"))
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
    })
}

#----- Main App UI & Server ------------------------------------------
ui <- fluidPage(
    titlePanel("Interactive WGCNA Explorer"),
    sidebarLayout(
        sidebarPanel(
            h4("Sample Clustering"),
            selectInput("distMethod", "Distance Method",
                choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                selected = "euclidean"
            ),
            sliderInput("cutHeight", "Cutoff Height:", min = 0, max = 50, value = 15, step = 1),
            hr(),
            h4("Scale-Free Topology"),
            sliderInput("powerRange", "Power Vector:", min = 1, max = 30, value = c(1, 20)),
            numericInput("rsqCut", "RsquaredCut:", value = 0.8, min = 0, max = 1, step = 0.05),
            numericInput("selectedPower", "Selected Power:", value = 6, min = 1, max = 30),
            hr(),
            h4("Module Detection"),
            numericInput("deepSplit", "deepSplit:", value = 2, min = 0, max = 4, step = 1),
            numericInput("minModuleSize", "Min Module Size:", value = 30, min = 5, max = 200, step = 1),
            actionButton("runModules", "Detect Modules"),
            width = 3
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("Sample Tree", sampleClustUI("sample")),
                tabPanel("Scale-Free Topology", sftUI("sft")),
                tabPanel("Gene Modules", geneModuleUI("mod")),
                tabPanel("Gene List", geneListUI("list"))
            ),
            width = 9
        )
    )
)

server <- function(input, output, session) {
  # Ensure expression.data is loaded
  validate(
    need(exists("expression.data"), "Load 'expression.data' before running app.")
  )

  # 1) Data quality check: goodSamplesGenes
  gsg <- goodSamplesGenes(expression.data, verbose = 0)
  if (!gsg$allOK) {
    # Filter out problematic samples/genes
    expression.data <- expression.data[gsg$goodSamples, gsg$goodGenes]
  }
  
  # 2) Numeric conversion with dimension validation
  exprDataNumeric <- reactive({
    mat <- expression.data
    df <- as.data.frame(lapply(mat, as.numeric))
    validate(
      need(nrow(df) > 1 && ncol(df) > 1, 
           "Filtered data has too few samples or genes to proceed.")
    )
    df
  })

  # 3) Call Shiny modules with cleaned data
  sampleClustServer(
    "sample", 
    exprData = exprDataNumeric, 
    distMethod = reactive(input$distMethod), 
    cutHeight = reactive(input$cutHeight)
  )

  sftServer(
    "sft", 
    exprData = exprDataNumeric, 
    powerRange = reactive(input$powerRange), 
    rsqCut = reactive(input$rsqCut), 
    selectedPower = reactive(input$selectedPower)
  )

  modulesObj <- geneModuleServer(
    "mod", 
    exprData = exprDataNumeric, 
    power = reactive(input$selectedPower), 
    deepSplit = reactive(input$deepSplit), 
    minSize = reactive(input$minModuleSize), 
    runTrigger = reactive(input$runModules)
  )

  geneListServer(
    "list", 
    exprData = exprDataNumeric, 
    modulesObj = modulesObj
  )
}

shinyApp(ui, server)
