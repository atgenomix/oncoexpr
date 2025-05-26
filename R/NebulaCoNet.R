#' @title shiny app for RNAseq WGCNA analysis
#' @description Simplified Shiny application keeping only data query (normcount, colData, DEG) and WGCNA modules
#' @import sparklyr
#' @import DBI
#' @import SummarizedExperiment
#' @import MultiAssayExperiment
#' @import shiny
#' @import shinyjs
#' @import shinybusy
#' @import DT
#' @import InteractiveComplexHeatmap
#' @import WGCNA
#' @importFrom DT dataTableOutput renderDataTable
#' @export
NebulaCoNet <- function(master = "sc://172.18.0.1:15002", method = "spark_connect", version = "3.5") {
  # Parallel plan
  future::plan(future::multisession, workers = 5)

  # UI definition
  ui <- fluidPage(
    navbarPage(
      title = "WGCNA App",
      # Data query tab
      tabPanel(
        title = "Gene Expression Profile",
        layout_sidebar(
          full_screen = TRUE,
          sidebar = sidebar(
            style = "min-height: 600px; overflow-y: auto;",
            h4("Analysis Runs"),
            dbBrowserUI("dbBrowser1")
          ),
          mainPanel(
            width = 12,
            tabsetPanel(
              tabPanel(
                "normCount Table",
                withSpinner(DT::dataTableOutput("wide_table_dt", width = "100%"))
              ),
              tabPanel(
                "DEG Table",
                downloadButton("download_DEG", "Download"),
                withSpinner(DT::dataTableOutput("DEG_table", width = "100%"))
              )
            ),
            progressPopupUI("popupProgress")
          )
        )
      ),
      # WGCNA tab
      tabPanel(
        title = "WGCNA",
        sidebarLayout(
          sidebarPanel(
            h4("Parallel Settings"),
            numericInput("nThreads", "Number of Threads:", value = parallel::detectCores() - 1, min = 1, max = parallel::detectCores(), step = 1),
            h4("Sample Clustering"),
            selectInput("distMethod", "Distance Method",
                        choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                        selected = "euclidean"),
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
            width = 9,
            tabsetPanel(
              tabPanel("Sample Tree", sampleClustUI("sample")),
              tabPanel("Scale-Free Topology", sftUI("sft")),
              tabPanel("Gene Modules", geneModuleUI("mod")),
              tabPanel("Gene List", geneListUI("list"))
            )
          )
        )
      )
    )
  )

  # Server logic
  server <- function(input, output, session) {
    # Reactive storage
    results <- reactiveValues(
      db_info = NULL,
      table_list = NULL,
      normcount_data = NULL,
      exacttest_data = NULL,
      coldata = NULL
    )

    # Spark connection
    sc <- sparklyr::spark_connect(master = master, method = method, version = version)
    session$onSessionEnded(function() {
      sparklyr::spark_disconnect(sc)
    })

    # Progress popup
    progressMod <- progressPopupServer("popupProgress")

    # Initialize DB browser and load table names
    observeEvent(sc, {
      results$db_info <- dbBrowserServer("dbBrowser1", sc)
    }, ignoreInit = FALSE)

    # When a database is selected: query tables, select latest normcount/exacttest/coldata
    observeEvent(results$db_info$selected_db(), {
      req(results$db_info$selected_db())
      dbName <- results$db_info$selected_db()
      # Stage 1: list tables and filter
      DBI::dbExecute(sc, paste0("USE ", dbName))
      tbls <- DBI::dbGetQuery(sc, paste0("SHOW TABLES IN ", dbName))$tableName
      prefixes <- c("^normcounts", "^exacttest", "^coldata")
      filtered <- tbls[grepl(paste(prefixes, collapse="|"), tbls, ignore.case=TRUE)]
      latest_df <- get_latest_file_group_df(filtered)
      sel <- latest_df$is_latest
      final_tbls <- latest_df$file[sel]
      norm_tbl <- final_tbls[grepl("^normcounts", final_tbls, ignore.case=TRUE)][1]
      exact_tbl <- final_tbls[grepl("^exacttest", final_tbls, ignore.case=TRUE)][1]
      col_tbl   <- final_tbls[grepl("^coldata",   final_tbls, ignore.case=TRUE)][1]

      # Stage 2: launch queries asynchronously
      norm_p <- future::future({
        conn <- sparklyr::spark_connect(master=master, method=method, version=version)
        on.exit(sparklyr::spark_disconnect(conn))
        DBI::dbExecute(conn, paste0("USE ", dbName))
        df <- DBI::dbGetQuery(conn, paste0("SELECT * FROM ", norm_tbl))
        names(df)[names(df)=="genes"] <- "GeneSymbol"
        df[ , names(df)!="_c0"]
      })
      exact_p <- future::future({
        conn <- sparklyr::spark_connect(master=master, method=method, version=version)
        on.exit(sparklyr::spark_disconnect(conn))
        DBI::dbExecute(conn, paste0("USE ", dbName))
        df <- DBI::dbGetQuery(conn, paste0("SELECT * FROM ", exact_tbl))
        names(df)[names(df)=="genes"] <- "GeneSymbol"
        df[ , names(df)!="_c0"]
      })
      col_p <- future::future({
        conn <- sparklyr::spark_connect(master=master, method=method, version=version)
        on.exit(sparklyr::spark_disconnect(conn))
        DBI::dbExecute(conn, paste0("USE ", dbName))
        DBI::dbGetQuery(conn, paste0("SELECT * FROM ", col_tbl))
      })

      progressMod$addPromise(norm_p, label="normcount")
      progressMod$addPromise(exact_p, label="exacttest")
      progressMod$addPromise(col_p,   label="coldata")

      # Stage 3: collect data
      promises::promise_all(
        norm = norm_p,
        ex = exact_p,
        col = col_p
      ) %...>% (function(res) {
        results$normcount_data <- res$norm
        results$exacttest_data <- res$ex
        results$coldata       <- res$col
      })
    })

    # Render normCount table
    output$wide_table_dt <- DT::renderDataTable({
      df <- results$normcount_data
      req(df)
      df_round <- as.data.frame(lapply(df, function(x) if(is.numeric(x)) round(x,4) else x))
      DT::datatable(df_round, options=list(pageLength=20, autoWidth=TRUE))
    })

    # Render DEG table and download handler
    output$DEG_table <- DT::renderDataTable({
      df <- results$exacttest_data
      req(df)
      df$logFC   <- round(df$logFC,5)
      df$logCPM  <- round(df$logCPM,5)
      df$PValue  <- formatC(df$PValue, format="e", digits=5)
      df$FDR     <- formatC(df$FDR,    format="e", digits=5)
      DT::datatable(df, options=list(pageLength=20, autoWidth=TRUE), server=FALSE)
    })
    output$download_DEG <- downloadHandler(
      filename = function() paste0("DEG_table_", Sys.Date(), ".csv"),
      content  = function(file) write.csv(results$exacttest_data, file, row.names=FALSE)
    )

    # Assemble SummarizedExperiment and MultiAssayExperiment
    observe({
      req(results$normcount_data, results$exacttest_data, results$coldata)
      assay_mat <- as.matrix(results$normcount_data[ , setdiff(names(results$normcount_data), "GeneSymbol")])
      rownames(assay_mat) <- results$normcount_data$GeneSymbol
      deg_mat   <- results$exacttest_data
      rownames(deg_mat) <- deg_mat$GeneSymbol
      col_tbl <- results$coldata
      rownames(col_tbl) <- names(assay_mat)
      se <- SummarizedExperiment(assays=list(normCount=assay_mat), colData=col_tbl, rowData=S4Vectors::DataFrame(deg_mat))
      mae <- MultiAssayExperiment(experiments=list(RNAseq=se), colData=col_tbl)
      # Store for WGCNA modules
      settingMAE <- reactiveVal(mae)
    })

    # Enable WGCNA threads
    observe({ req(input$nThreads); enableWGCNAThreads(input$nThreads) })

    # WGCNA analysis pipeline
    observeEvent(results$normcount_data, {
      df <- results$normcount_data
      req(df)
      rownames(df) <- df$GeneSymbol
      expr <- t(as.matrix(df[ , setdiff(names(df), "GeneSymbol")]))
      # Quality check
      gsg <- WGCNA::goodSamplesGenes(expr, verbose=0)
      if(!gsg$allOK) expr <- expr[gsg$goodSamples, gsg$goodGenes]
      # Reactive expression for cleaned data
      exprDataNumeric <- reactive({ as.data.frame(lapply(as.data.frame(expr), as.numeric)) })
      # Call WGCNA modules
      sampleOut <- sampleClustServer("sample", exprData=exprDataNumeric, distMethod=reactive(input$distMethod), cutHeight=reactive(input$cutHeight))
      observeEvent(sampleOut$maxHeight(), {
        updateSliderInput(session, "cutHeight", min=0, max=ceiling(sampleOut$maxHeight()), value=min(input$cutHeight, ceiling(sampleOut$maxHeight())))
      })
      sftServer("sft", exprData=sampleOut$filteredExpr, powerRange=reactive(input$powerRange), rsqCut=reactive(input$rsqCut), selectedPower=reactive(input$selectedPower))
      modulesObj <- geneModuleServer("mod", exprData=sampleOut$filteredExpr, power=reactive(input$selectedPower), deepSplit=reactive(input$deepSplit), minSize=reactive(input$minModuleSize), runTrigger=reactive(input$runModules))
      geneListServer("list", exprData=sampleOut$filteredExpr, modulesObj=modulesObj)
    })
  }

  # Run the app
  shinyApp(ui=ui, server=server)
}
