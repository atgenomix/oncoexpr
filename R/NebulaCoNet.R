#' @import edgeR
#' @import limma
#' @import MultiAssayExperiment
#' @import SummarizedExperiment
#' @import pcaPP
#' @import reshape2
#' @import stringr
#' @import readxl
#' @import ggplot2
#' @import tidyr
#' @import tibble
#' @import viridis
#' @import RColorBrewer
#' @import pheatmap
#' @import ggrepel
#' @import readr
#' @import dplyr
#' @import sparklyr
#' @import org.Hs.eg.db
#' @import enrichplot
#' @import clusterProfiler
#' @import DBI
#' @import shinybusy
#' @import shiny
#' @import bslib
#' @import DT
#' @import InteractiveComplexHeatmap
#' @import ggiraph
#' @import shinydashboard
#' @import promises
#' @import future
#' @import purrr
#' @import shinycssloaders
#' @import shinyjs
#' @import factoextra
#' @import ggforce
#' @import ggrepel
#' @importFrom ggpubr color_palette
#' @importFrom enrichplot color_palette
#' @importFrom DT dataTableOutput renderDataTable

#'
NULL
#' @keywords internal
#' RNAseqShinyAppSpark
#'
#' This function creates a Shiny application for RNA sequencing data analysis using Spark.
#' It provides an interactive interface for users to explore and analyze RNAseq data.
#'

#' @title shiny app for RNAseq for public use
#' @description start the RNAseq shiny app with spark connection for shiny proxy at the seqslab console
#' @return A Shiny application object.
#'
#' @examples
#' \dontrun{
#' library(sparklyr)
#' sc <- spark_connect(master = "local")
#' input_data <- read.csv("path/to/rnaseq_data.csv")
#' RNAseqShinyAppSpark(input_data, sc, "RNAseq Analysis", 8080)
#' }
#'
#' @export


NebulaCoNet <- function(master = "sc://172.18.0.1:15002", method = "spark_connect", version = "3.5") {
  plan(multisession, workers = 5)
  # plan(sequential)
  print(future::plan())
  ui <- fluidPage(
    navbarPage(
      title = "NebulaCoNet",
      tabPanel(
        title = "Data Query",
        layout_sidebar(
          full_screen = TRUE,
          sidebar = sidebar(
            style = "min-height: 600px; overflow-y: auto;",
            h4("Analysis Runs"),
            dbBrowserUI("dbBrowser1")
          ),
          mainPanel(
            fluidRow(
              useShinyjs(),
              column(
                width = 12,
                tabsetPanel(
                  tabPanel(
                    "normCount Table",
                    withSpinner(DT::dataTableOutput("wide_table_dt", width = "100%"))
                  ),
                  tabPanel(
                    "DEG Table",
                    downloadButton("download_DEG", "Download"),
                    withSpinner((DT::dataTableOutput("DEG_table", width = "100%")))
                  )
                )
              ),
              progressPopupUI("popupProgress")
            )
          )
        )
      ),
      tabPanel(
        title = "WGCNA",
        sidebarLayout(
          sidebarPanel(
              h4("Parallel Settings"),
              numericInput(
                "nThreads", "Number of Threads:",
                value = parallel::detectCores() - 1,
                min = 1,
                max = parallel::detectCores(),
                step = 1
              ),
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
    )
  )

  server <- function(input, output, session) {
    results <- reactiveValues(
      db_info = NULL,
      table_list = NULL,
      grouplist = NULL,
      normcount_data = NULL,
      exacttest_data = NULL,
      coldata = NULL
    )
    volcano_res <- reactiveVal(NULL)
    settingMAE <- reactiveVal(NULL)
    DEG_table <- reactiveVal(NULL)
    DEG_summary <- reactiveVal(NULL)
    wide_data <- reactiveVal(NULL)
    maeColData <- reactiveVal(NULL)


    sc <- sparklyr::spark_connect(master = master, method = method, version = version)

    session$onSessionEnded(function() {
      if (!is.null(sc)) {
        sparklyr::spark_disconnect(sc)
        message("Spark connection disconnected.")
      }
    })

    observeEvent(sc,
      {
        req(sc)
        print("dbbrowser initialized")
        shinyjs::disable("dbBrowser1-selected_db")
        results$db_info <- dbBrowserServer("dbBrowser1", sc)
        showNotification("Waiting for initialization", type = "message", duration = 10)

        normcount_future <- trigger_cluster_query_by_pattern(
          master, method, version,
          pattern = "^normcounts",
          output_label = "init_tbl_normcount"
        )

        exacttest_future <- trigger_cluster_query_by_pattern(
          master, method, version,
          pattern = "^exacttest",
          output_label = "init_tbl_exacttest"
        )

        coldata_future <- trigger_cluster_query_by_pattern(
          master, method, version,
          pattern = "^coldata",
          output_label = "init_tbl_coldata"
        )
        all_promises <- promises::promise_all(
          norm = normcount_future,
          ex   = exacttest_future,
          col  = coldata_future
        )

        all_promises %...>% (function(res_list) {
          shinyjs::enable("dbBrowser1-selected_db")
          showNotification("Initialization complete. Check list!", type = "message", duration = 10)
        }) %...!% (function(e) {
          shinyjs::enable("dbBrowser1-selected_db")
          showNotification(paste("Error:", e$message), type = "error")
        })

      },
      ignoreInit = FALSE
    )

    progressMod <- progressPopupServer("popupProgress")

    observeEvent(results$db_info$selected_db(), {
      req(results$db_info$selected_db())
      shinyjs::disable("dbBrowser1-selected_db")
      output$wide_table_dt <- DT::renderDataTable({
        data.frame()
      })
      output$DEG_table <- DT::renderDataTable({
        data.frame()
      })

      results$exacttest_data <- NULL
      results$normcount_data <- NULL
      results$coldata <- NULL
      output$ht_heatmap <- renderPlot({
        grid::grid.newpage()
        grid::grid.text("No data available.")
      })

      selected_db_name <- results$db_info$selected_db()
      message(sprintf("[DB Selected] %s at %s", selected_db_name, Sys.time()))

      withProgress(message = "Stage 1: Listing & filtering tables", value = 0, {
        t0 <- Sys.time()
        message(sprintf("[Stage1] Start at %s", t0))

        DBI::dbExecute(sc, paste0("USE ", selected_db_name))
        tbl_list_query <- DBI::dbGetQuery(sc, paste0("SHOW TABLES IN ", selected_db_name))
        tbls <- tbl_list_query$tableName
        t1 <- Sys.time()
        message(sprintf(
          "[Stage1] Fetched %d tables at %s (%.2f sec)",
          length(tbls), t1, as.numeric(difftime(t1, t0, "secs"))
        ))
        setProgress(value = 0.2, detail = sprintf("Fetched %d tables", length(tbls)))

        prefix <- c("^normcounts|^exacttest|^coldata")
        tbl_list_query_prefix <- tbl_list_query[grepl(paste(prefix, collapse = "|"), tbls), ]
        t2 <- Sys.time()
        message(sprintf(
          "[Stage1] Prefix filter → %d tables at %s (%.2f sec)",
          nrow(tbl_list_query_prefix), t2, as.numeric(difftime(t2, t1, "secs"))
        ))
        setProgress(value = 0.4, detail = sprintf("Prefix filter: %d tables", nrow(tbl_list_query_prefix)))

        tbls_with_prefix <- tbl_list_query_prefix$tableName
        tbls_with_time_filter <- get_latest_file_group_df(tbls_with_prefix)
        t3 <- Sys.time()
        message(sprintf(
          "[Stage1] Time filter applied at %s (%.2f sec)",
          t3, as.numeric(difftime(t3, t2, "secs"))
        ))
        setProgress(value = 0.6, detail = "Applied time filter")

        if (any(tbls_with_time_filter$is_latest)) {
          sel_idx <- tbls_with_time_filter$is_latest
          message(sprintf("[Stage1] Latest tables found at %s", Sys.time()))
        } else {
          sel_idx <- !tbls_with_time_filter$is_latest
          message(sprintf("[Stage1] No latest table, using all at %s", Sys.time()))
        }
        tbl_list_query_prefix_time <- tbl_list_query_prefix[sel_idx, ]
        summary_table <- tbls_with_time_filter[sel_idx, ]
        t4 <- Sys.time()
        message(sprintf(
          "[Stage1] Selected %d tables at %s (%.2f sec)",
          nrow(tbl_list_query_prefix_time), t4, as.numeric(difftime(t4, t3, "secs"))
        ))
        setProgress(value = 0.8, detail = sprintf("Selected %d tables", nrow(tbl_list_query_prefix_time)))

        tbls_final <- summary_table$file
        normcount_tbls <- tbl_list_query_prefix_time[grepl("^normcounts", tbls_final, ignore.case = TRUE), "tableName"]
        exacttest_tbls <- tbl_list_query_prefix_time[grepl("^exacttest", tbls_final, ignore.case = TRUE), "tableName"]
        coldata_tbls <- tbl_list_query_prefix_time[grepl("^coldata", tbls_final, ignore.case = TRUE), "tableName"]
        t5 <- Sys.time()
        message(sprintf(
          "[Stage1] Categorized tables at %s (%.2f sec)",
          t5, as.numeric(difftime(t5, t4, "secs"))
        ))

        setProgress(value = 1, detail = "Stage 1 complete")
        results$table_list <- tbl_list_query_prefix_time
      })



      req(normcount_tbls, exacttest_tbls, coldata_tbls)

      t0_norm_launch <- Sys.time()
      message(sprintf("[Stage2-normcount] Launch at %s", t0_norm_launch))
      normcount_promise <- future_promise(
        {
          start_time <- Sys.time()
          message(sprintf("[%s] Start querying normcounts table", start_time))

          sc_conn <- sparklyr::spark_connect(master = master, method = method, version = version)
          on.exit(sparklyr::spark_disconnect(sc_conn))
          DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
          query_normcount <- paste0("SELECT * FROM ", normcount_tbls[1])
          normcount <- DBI::dbGetQuery(sc_conn, query_normcount)

          colnames(normcount)[colnames(normcount) == "genes"] <- "GeneSymbol"
          normcount <- normcount[, colnames(normcount) != "_c0"]

          end_time <- Sys.time()
          message(sprintf(
            "[%s] Completed normcounts query (Duration: %.2f seconds)",
            end_time, as.numeric(difftime(end_time, start_time, units = "secs"))
          ))

          normcount
        },
        globals = list(
          master = master, method = method, version = version,
          normcount_tbls = normcount_tbls, selected_db_name = selected_db_name
        ),
        seed = TRUE
      )

      t0_exact_launch <- Sys.time()
      message(sprintf("[Stage2-exacttest] Launch at %s", t0_exact_launch))
      exacttest_promise <- future_promise(
        {
          start_time <- Sys.time()
          message(sprintf("[%s] Start querying exacttest table", start_time))

          sc_conn <- sparklyr::spark_connect(master = master, method = method, version = version)
          on.exit(sparklyr::spark_disconnect(sc_conn))
          DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
          query_exacttest <- paste0("SELECT * FROM ", exacttest_tbls[1])
          exacttest <- DBI::dbGetQuery(sc_conn, query_exacttest)

          colnames(exacttest)[colnames(exacttest) == "genes"] <- "GeneSymbol"
          exacttest <- exacttest[, colnames(exacttest) != "_c0"]

          end_time <- Sys.time()
          message(sprintf(
            "[%s] Completed exacttest query (Duration: %.2f seconds)",
            end_time, as.numeric(difftime(end_time, start_time, units = "secs"))
          ))

          exacttest
        },
        globals = list(
          master = master, method = method, version = version,
          exacttest_tbls = exacttest_tbls, selected_db_name = selected_db_name
        ),
        seed = TRUE
      )

      t0_coldata_launch <- Sys.time()
      message(sprintf("[Stage2-coldata] Launch at %s", t0_coldata_launch))
      coldata_promise <- if (length(coldata_tbls) > 0) {
        future_promise(
          {
            start_time <- Sys.time()
            message(sprintf("[%s] Start querying coldata table", start_time))

            sc_conn <- sparklyr::spark_connect(master = master, method = method, version = version)
            on.exit(sparklyr::spark_disconnect(sc_conn))
            DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
            query_coldata <- paste0("SELECT * FROM ", coldata_tbls[1])
            coldata <- DBI::dbGetQuery(sc_conn, query_coldata)

            end_time <- Sys.time()
            message(sprintf(
              "[%s] Completed coldata query (Duration: %.2f seconds)",
              end_time, as.numeric(difftime(end_time, start_time, units = "secs"))
            ))

            coldata
          },
          globals = list(
            master = master, method = method, version = version,
            coldata_tbls = coldata_tbls, selected_db_name = selected_db_name
          ),
          seed = TRUE
        )
      } else {
        normcount_promise %...>% (function(normcount) {
          future_promise(
            {
              start_time <- Sys.time()
              message(sprintf("[%s] Generating random coldata", start_time))
              coldata <- generate_colData_random(normcount, genecol = "GeneSymbol")
              end_time <- Sys.time()
              message(sprintf(
                "[%s] Completed coldata generation (Duration: %.2f seconds)",
                end_time, as.numeric(difftime(end_time, start_time, units = "secs"))
              ))
              coldata
            },
            seed = TRUE
          )
        })
      }
      progressMod$addPromise(normcount_promise, label = "normcount")
      progressMod$addPromise(coldata_promise, label = "coldata")
      progressMod$addPromise(exacttest_promise, label = "exacttest")


      withProgress(message = "Stage 3: Collecting and processing data", value = 0, {
        t0_collect <- Sys.time()
        message(sprintf("[Stage3] Collection start at %s", t0_collect))

        promise_all(
          normcount_data = normcount_promise,
          exacttest_data = exacttest_promise,
          coldata        = coldata_promise
        ) %...>% with({
          t1_collect <- Sys.time()
          message(sprintf(
            "[Stage3] Collection completed at %s (Duration: %.2f seconds)",
            t1_collect, as.numeric(difftime(t1_collect, t0_collect, units = "secs"))
          ))

          results$normcount_data <- normcount_data
          results$exacttest_data <- exacttest_data
          results$coldata <- coldata

          message("=== normcount_data ===")
          print(head(results$normcount_data))
          message("=== exacttest_data ===")
          print(head(results$exacttest_data))
          message("=== coldata ===")
          print(head(results$coldata))

          setProgress(value = 1, detail = "Data collected")
        })
      })

      wide_data(NULL)
      DEG_table(NULL)

      observe({
        req(results$exacttest_data, results$normcount_data, results$coldata)
        DEG_table(results$exacttest_data)
        wide_data(results$normcount_data)
        maeColData(results$coldata)
        message("assign reactiveVal: DEG_table, wide_data, maeColData")
      })
      output$wide_table_dt <- DT::renderDataTable({
        req(wide_data())
        normCount_round <- as.data.frame(lapply(
          wide_data(),
          function(x) if (is.numeric(x)) round(x, 4) else x
        ))

        print("send wide data to UI")
        DT::datatable(
          normCount_round,
          options = list(pageLength = 20, autoWidth = TRUE)
        )
      })

      message(sprintf("[Process] Completed all stages on PID %s at %s", Sys.getpid(), Sys.time()))
    })

    observe({
      deg_table_round <- DEG_table()
      req(deg_table_round)
      deg_table_round$logFC <- if (is.numeric(deg_table_round$logFC)) round(deg_table_round$logFC, 5) else deg_table_round$logFC
      deg_table_round$logCPM <- if (is.numeric(deg_table_round$logCPM)) round(deg_table_round$logCPM, 5) else deg_table_round$logCPM
      deg_table_round$PValue <- formatC(deg_table_round$PValue, format="e", digits = 5)
      deg_table_round$FDR <- formatC(deg_table_round$FDR, format="e", digits = 5)

      output$DEG_table <- renderDT(
        {
          datatable(
            deg_table_round,
            options = list(pageLength = 20, autoWidth = TRUE)
          )
        },
        server = FALSE
      )
    })

    observe({
      req(input$nThreads)
      enableWGCNAThreads(input$nThreads)
    })

    observeEvent(wide_data(), {
      req(wide_data())
        enableWGCNAThreads()
        raw_df <- wide_data()
        if (!"GeneSymbol" %in% colnames(raw_df)) stop("wide_data() must contain 'GeneSymbol' column")
        rownames(raw_df) <- raw_df$GeneSymbol
        expr_mat_raw <- raw_df[, setdiff(colnames(raw_df), "GeneSymbol"), drop = FALSE]
        print(dim(expr_mat_raw))
        # Filter numeric columns only
        numeric_df <- raw_df[, sapply(raw_df, is.numeric), drop = FALSE]
        # Transpose to samples x genes matrix
        expr_mat <- t(as.matrix(expr_mat_raw))
        # option for development
        #expr_mat <- as.data.frame(expr_mat[, sample(1:ncol(expr_mat), 2500, replace = FALSE)])
        print(dim(expr_mat))
        # 1) Data quality check: require numeric matrix
        gsg <- goodSamplesGenes(expr_mat, verbose = 0)
        print(gsg$allOK)
        if (!gsg$allOK) {
          expr_mat <- expr_mat[gsg$goodSamples, gsg$goodGenes]
        }
        
        # 2) Numeric conversion with dimension validation
        exprDataNumeric <- reactive({
          mat <- as.data.frame(expr_mat)
          df <- as.data.frame(lapply(mat, as.numeric))
          validate(
            need(nrow(df) > 1 && ncol(df) > 1, 
                "Filtered data has too few samples or genes to proceed.")
          )
          df
        })

        # 3) Call Shiny modules with cleaned data
        sampleOut <- sampleClustServer(
          "sample", 
          exprData = exprDataNumeric, 
          distMethod = reactive(input$distMethod), 
          cutHeight = reactive(input$cutHeight)
        )
        
        observeEvent(sampleOut$maxHeight(), {
          maxH <- sampleOut$maxHeight()
          updateSliderInput(
            session,
            "cutHeight",
            min = 0,
            max = ceiling(maxH),
            step = 1,
            value = min(input$cutHeight, ceiling(maxH))
          )
        })

        sftServer(
          "sft", 
          exprData = sampleOut$filteredExpr, 
          powerRange = reactive(input$powerRange), 
          rsqCut = reactive(input$rsqCut), 
          selectedPower = reactive(input$selectedPower)
        )

        modulesObj <- geneModuleServer(
          "mod", 
          exprData = sampleOut$filteredExpr, 
          power = reactive(input$selectedPower), 
          deepSplit = reactive(input$deepSplit), 
          minSize = reactive(input$minModuleSize), 
          runTrigger = reactive(input$runModules)
        )

        geneListServer(
          "list", 
          exprData = sampleOut$filteredExpr, 
          modulesObj = modulesObj
        )
    })



  }

  for_run <- shinyApp(ui = ui, server = server)
  runapp <- runApp(for_run)

  return(runapp)
}
