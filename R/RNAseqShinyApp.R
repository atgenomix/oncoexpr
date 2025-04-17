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


RNAseqShinyAppSpark <- function(master = "sc://172.18.0.1:15002", method = "spark_connect", version = "3.5") {
  plan(multisession, workers = 4)
  # plan(sequential)
  print(future::plan())
  ui <- fluidPage(
    navbarPage(
      title = "RNAseq App",
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
                    downloadButton("download_DEG", "Download DEG CSV"),
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
        title = "Principal Component Analysis",
        mainPanel(
          width = 12,
          pcaModuleUI("pca1")
        )
      ),
      tabPanel(
        title = "Differential Expression Analysis",
        layout_sidebar(
          full_screen = TRUE,
          sidebar = sidebar(
            style = "min-height: 600px; overflow-y: auto;",
            sliderInput("lfc_cut", "Fold Change Threshold (log2):",
              min = 0, max = 10, value = 1, step = 0.1
            ),
            sliderInput("pval_cut", "p-value Threshold:",
              min = 0.001, max = 1, value = 0.05, step = 0.001
            ),
            sliderInput("pointSize", "Point Size:",
              min = 1, max = 5, value = 2, step = 0.5
            ),
            sliderInput("ptAlpha", "Transparent:",
              min = 0.1, max = 1, value = 0.6, step = 0.1
            ),
            sliderInput("labelSize", "Gene Label Size:",
              min = 1, max = 6, value = 3, step = 0.5
            ),
            numericInput("topN", "Label N number of Genes (0 is no labeling):",
              value = 100, min = 0, max = 1000
            ),
            actionButton("run_DEG", "Update DEG")
          ),
          mainPanel(
            width = 12,
            tabsetPanel(
              tabPanel(
                "Volcano Plot",
                fluidRow(
                  column(
                    width = 12,
                    interactivePlotsUI("volcano_plots")
                  ),
                  column(
                    width = 12,
                    mod_geneSelector_ui("gene_selector")
                  )
                )
              ),
              tabPanel(
                "Heatmap",
                fluidRow(
                  column(
                    width = 6,
                    box(
                      title = "Differential heatmap", width = NULL, solidHeader = TRUE, status = "primary",
                      originalHeatmapOutput("ht", height = 1000, containment = TRUE)
                    )
                  ),
                  column(
                    width = 6,
                    box(
                      title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
                      subHeatmapOutput("ht", title = NULL, containment = TRUE)
                    )
                  )
                )
              ),
              tabPanel(
                "Gene Set Enrichment",
                fluidRow(
                  column(
                    12,
                    h4("UPregulated DEGs"),
                    tabsetPanel(
                      tabPanel("KEGG", withSpinner(plotOutput("G1_KEGG"))),
                      tabPanel("MF", withSpinner(plotOutput("G1_MF"))),
                      tabPanel("BP", withSpinner(plotOutput("G1_BP"))),
                      tabPanel("CC", withSpinner(plotOutput("G1_CC"))),
                      tabPanel("GSEA(KEGG)", withSpinner(gseaFCModuleUI("gsea_up")))
                    )
                  )
                ),
                fluidRow(
                  column(
                    12,
                    h4("DOWNregulated DEGs"),
                    tabsetPanel(
                      tabPanel("KEGG", withSpinner(plotOutput("G2_KEGG"))),
                      tabPanel("MF", withSpinner(plotOutput("G2_MF"))),
                      tabPanel("BP", withSpinner(plotOutput("G2_BP"))),
                      tabPanel("CC", withSpinner(plotOutput("G2_CC"))),
                      tabPanel("GSEA(KEGG)", withSpinner(gseaFCModuleUI("gsea_down")))
                    )
                  )
                )
              )
            )
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

    # init executors 
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

    #query tables from the database
    observeEvent(results$db_info$selected_db(), {
      req(results$db_info$selected_db())
      shinyjs::disable("dbBrowser1-selected_db")

      selected_db_name <- results$db_info$selected_db()
      message(sprintf("[Stage0] DB selected → %s at %s",
                      selected_db_name, Sys.time()))

      # Stage 1: Listing & filtering tables
      t_stage1 <- system.time({
        DBI::dbExecute(sc, paste0("USE ", selected_db_name))
        tbl_list_query <- DBI::dbGetQuery(sc,
          paste0("SHOW TABLES IN ", selected_db_name))
        tbls <- tbl_list_query$tableName

        prefix <- c("^normcounts", "^exacttest", "^coldata")
        tbl_pref    <- tbl_list_query[grepl(paste(prefix, collapse="|"), tbls), ]
        tbls_time   <- get_latest_file_group_df(tbl_pref$tableName)
        sel         <- if (any(tbls_time$is_latest)) tbls_time$is_latest else rep(TRUE, nrow(tbls_time))
        tbl_final_meta <- tbl_pref[sel, ]
      })
      message(sprintf("[Stage1] Listing & filtering complete at %s (%.2f s)",
                      Sys.time(), t_stage1["elapsed"]))

      normcount_tbls <- tbl_final_meta$tableName[
        grepl("^normcounts", tbl_final_meta$tableName, ignore.case=TRUE)]
      exacttest_tbls <- tbl_final_meta$tableName[
        grepl("^exacttest", tbl_final_meta$tableName, ignore.case=TRUE)]
      coldata_tbls   <- tbl_final_meta$tableName[
        grepl("^coldata", tbl_final_meta$tableName, ignore.case=TRUE)]

      # Stage 2: Launching queries
      t2_norm <- system.time({
        normcount_promise <- future_promise({
          start_time <- Sys.time()
          message(sprintf("[%s] Start querying normcounts table", start_time))

          sc_conn <- sparklyr::spark_connect(master = master,
                                              method = method,
                                              version = version)
          on.exit(sparklyr::spark_disconnect(sc_conn))
          DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
          query_norm <- paste0("SELECT * FROM ", normcount_tbls[1])
          normcount <- DBI::dbGetQuery(sc_conn, query_norm)
          colnames(normcount)[colnames(normcount)=="genes"] <- "GeneSymbol"
          normcount <- normcount[, colnames(normcount)!="_c0"]

          end_time <- Sys.time()
          message(sprintf("[%s] Completed normcounts (%.2f s)",
                          end_time,
                          as.numeric(difftime(end_time, start_time, "secs"))))
          normcount
        }, globals = list(
          master = master, method = method, version = version,
          normcount_tbls = normcount_tbls,
          selected_db_name = selected_db_name
        ), seed = TRUE)
      })
      progressMod$addPromise(normcount_promise, label = "normcount")
      message(sprintf("[Stage2-norm] Launched at %s (%.2f s)",
                      Sys.time(), t2_norm["elapsed"]))

      t2_exact <- system.time({
        exacttest_promise <- future_promise({
          start_time <- Sys.time()
          message(sprintf("[%s] Start querying exacttest table", start_time))

          sc_conn <- sparklyr::spark_connect(master = master,
                                              method = method,
                                              version = version)
          on.exit(sparklyr::spark_disconnect(sc_conn))
          DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
          query_exact <- paste0("SELECT * FROM ", exacttest_tbls[1])
          exacttest <- DBI::dbGetQuery(sc_conn, query_exact)
          colnames(exacttest)[colnames(exacttest)=="genes"] <- "GeneSymbol"
          exacttest <- exacttest[, colnames(exacttest)!="_c0"]

          end_time <- Sys.time()
          message(sprintf("[%s] Completed exacttest (%.2f s)",
                          end_time,
                          as.numeric(difftime(end_time, start_time, "secs"))))
          exacttest
        }, globals = list(
          master = master, method = method, version = version,
          exacttest_tbls = exacttest_tbls,
          selected_db_name = selected_db_name
        ), seed = TRUE)
      })
      progressMod$addPromise(exacttest_promise, label = "exacttest")
      message(sprintf("[Stage2-exact] Launched at %s (%.2f s)",
                      Sys.time(), t2_exact["elapsed"]))

      t2_coldata <- system.time({
        coldata_promise <- if (length(coldata_tbls) > 0) {
          future_promise({
            start_time <- Sys.time()
            message(sprintf("[%s] Start querying coldata table", start_time))

            sc_conn <- sparklyr::spark_connect(master = master,
                                                method = method,
                                                version = version)
            on.exit(sparklyr::spark_disconnect(sc_conn))
            DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
            query_col <- paste0("SELECT * FROM ", coldata_tbls[1])
            coldata <- DBI::dbGetQuery(sc_conn, query_col)

            end_time <- Sys.time()
            message(sprintf("[%s] Completed coldata (%.2f s)",
                            end_time,
                            as.numeric(difftime(end_time, start_time, "secs"))))
            coldata
          }, globals = list(
            master = master, method = method, version = version,
            coldata_tbls = coldata_tbls,
            selected_db_name = selected_db_name
          ), seed = TRUE)
        } else {
          normcount_promise %...>% (function(nc) future_promise({
            start_time <- Sys.time()
            message(sprintf("[%s] Generating random coldata", start_time))
            cd <- generate_colData_random(nc, genecol="GeneSymbol")
            end_time <- Sys.time()
            message(sprintf("[%s] Completed random coldata (%.2f s)",
                            end_time,
                            as.numeric(difftime(end_time, start_time, "secs"))))
            cd
          }, seed = TRUE))
        }
      })
      progressMod$addPromise(coldata_promise, label = "coldata")
      message(sprintf("[Stage2-coldata] Launched at %s (%.2f s)",
                      Sys.time(), t2_coldata["elapsed"]))

      # Stage 3: Collecting data
      t_stage3 <- system.time({
        promises::promise_all(
          normcount_data = normcount_promise,
          exacttest_data = exacttest_promise,
          coldata        = coldata_promise
        ) %...>% with({
          results$normcount_data <- normcount_data
          results$exacttest_data <- exacttest_data
          results$coldata        <- coldata
        })
      })
      message(sprintf("[Stage3] Data collection complete at %s (%.2f s)",
                      Sys.time(), t_stage3["elapsed"]))

      # Stage 4: Assigning reactive values
      t_stage4 <- system.time({
        observe({
          req(results$exacttest_data, results$normcount_data, results$coldata)
          DEG_table(results$exacttest_data)
          wide_data(results$normcount_data)
          maeColData(results$coldata)
        })
      })
      message(sprintf("[Stage4] Reactive values assigned at %s (%.2f s)",
                      Sys.time(), t_stage4["elapsed"]))

      # Stage 5: UI update
      t_stage5 <- system.time({
        output$wide_table_dt <- DT::renderDataTable({
          req(wide_data())
          normCount_round <- as.data.frame(lapply(
            wide_data(),
            function(x) if (is.numeric(x)) round(x, 4) else x
          ))
          DT::datatable(normCount_round,
                        options = list(pageLength = 20, autoWidth = TRUE))
        })
        output$DEG_table <- DT::renderDataTable({
          req(DEG_table())
          deg_round <- DEG_table()
          deg_round$logFC  <- if (is.numeric(deg_round$logFC)) round(deg_round$logFC, 5) else deg_round$logFC
          deg_round$logCPM <- if (is.numeric(deg_round$logCPM)) round(deg_round$logCPM, 5) else deg_round$logCPM
          deg_round$PValue <- formatC(deg_round$PValue, format="e", digits=5)
          deg_round$FDR    <- formatC(deg_round$FDR,    format="e", digits=5)
          DT::datatable(deg_round,
                        options = list(pageLength = 20, autoWidth = TRUE))
        })
        shinyjs::enable("dbBrowser1-selected_db")
      })
      message(sprintf("[Stage5] UI update complete at %s (%.2f s)",
                      Sys.time(), t_stage5["elapsed"]))

      # Total
      total_elapsed <- t_stage1["elapsed"] +
                      t2_norm["elapsed"] + t2_exact["elapsed"] +
                      t2_coldata["elapsed"] + t_stage3["elapsed"] +
                      t_stage4["elapsed"] + t_stage5["elapsed"]
      message(sprintf("[Total] All stages done in %.2f seconds", total_elapsed))
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
      req(DEG_table(), wide_data(), maeColData())
      shinyjs::enable("dbBrowser1-selected_db")
      assay_data <- as.matrix(wide_data()[, -which(colnames(wide_data()) == "GeneSymbol")])
      if ("GeneSymbol" %in% colnames(wide_data())) {
        rownames(assay_data) <- wide_data()[, "GeneSymbol"]
      }

      deg_data <- DEG_table()
      if ("GeneSymbol" %in% colnames(deg_data)) {
        rownames(deg_data) <- deg_data$GeneSymbol
      }

      output$download_DEG <- downloadHandler(
        filename = function() {
          paste("DEG_table_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(DEG_table(), file, row.names = FALSE)
        }
      )
      common_genes <- intersect(rownames(assay_data), rownames(deg_data))
      assay_data <- assay_data[common_genes, , drop = FALSE]
      deg_data_sub <- deg_data[common_genes, , drop = FALSE]
      sample_info_table <- maeColData()
      rownames(sample_info_table) <- colnames(assay_data)

      se_expression_matrix <- SummarizedExperiment(
        assays = list(normCount = assay_data),
        colData = sample_info_table,
        rowData = S4Vectors::DataFrame(deg_data_sub)
      )

      mae <- MultiAssayExperiment(
        experiments = list(
          RNAseq = se_expression_matrix
        ),
        colData = sample_info_table
      )
      settingMAE(mae)
    })


    topGeneList <- reactiveVal(NULL)
    downGeneList <- reactiveVal(NULL)

    observe({
      req(DEG_table(), maeColData(), wide_data())
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


      normCount <- wide_data()
      volcanoData <- DEG_table()
      colData <- maeColData()

      exprData <- transfExprFormat(normCount, colData)


      DEG_table_data <- DEG_table()
      topGenes <- DEG_table_data[DEG_table_data$PValue < input$pval_cut & DEG_table_data$logFC > input$lfc_cut, "GeneSymbol"]
      downGenes <- DEG_table_data[DEG_table_data$PValue < input$pval_cut & DEG_table_data$logFC < -input$lfc_cut, "GeneSymbol"]
      topGeneList(topGenes)
      downGeneList(downGenes)

      geneListVec <- c(topGenes, downGenes)
      if (!is.null(geneListReactive)) {
        selected_gene <- mod_geneSelector_server("gene_selector", volcanoData, geneListVec)
      } else {
        (selected_gene <- NULL)
      }
      interactivePlotsServer("volcano_plots", volcanoData = volcanoData, exprData = exprData, params = params, selectedGene = selected_gene)
    })

    geneListReactive <- eventReactive(input$run_DEG, {
      req(DEG_table(), maeColData(), wide_data())
      shinyjs::disable("dbBrowser1-selected_db")
      output$ht_heatmap <- renderPlot({
        grid::grid.newpage()
        grid::grid.text("No data available.")
      })

      DEG_table_data <- DEG_table()
      topGenes <- DEG_table_data[DEG_table_data$PValue < input$pval_cut & DEG_table_data$logFC > input$lfc_cut, "GeneSymbol"]
      downGenes <- DEG_table_data[DEG_table_data$PValue < input$pval_cut & DEG_table_data$logFC < -input$lfc_cut, "GeneSymbol"]

      topGeneList(topGenes)
      downGeneList(downGenes)
      print("=====Case DEG List=====")
      print(length(topGeneList()))
      print("=====Control DEG List=====")
      print(length(downGeneList()))

      gene_list <- paste(c(topGenes, downGenes), collapse = ",")
      gene_list
    })

    observeEvent(input$run_DEG, {
      req(topGeneList(), downGeneList())
      shinyjs::disable("run_DEG")
    })

    observeEvent(geneListReactive(), {
      req(geneListReactive(), settingMAE())
      mae <- settingMAE()

      geneListVec <- unlist(strsplit(geneListReactive(), ","))
      geneListVec <- trimws(geneListVec)

      future_promise({
        make_heatmap_mae(mae, geneListVec)
      }) %...>% (function(ht) {
        if (!is.null(ht)) {
          makeInteractiveComplexHeatmap(input, output, session, ht, "ht")
        } else {
          output$ht_heatmap <- renderPlot({
            grid::grid.newpage()
            grid::grid.text("No data available.")
          })
        }
        outputOptions(output, "ht_heatmap", suspendWhenHidden = FALSE)
      }) %...!% (function(e) {
        output$ht_heatmap <- renderPlot({
          grid::grid.newpage()
          grid::grid.text(paste("An error occurred:", e$message))
        })
        outputOptions(output, "ht_heatmap", suspendWhenHidden = FALSE)
      })
    })

    result_G1_CC <- reactiveVal(NULL)
    result_G1_BP <- reactiveVal(NULL)
    result_G1_MF <- reactiveVal(NULL)
    result_G1_KEGG <- reactiveVal(NULL)

    result_G2_CC <- reactiveVal(NULL)
    result_G2_BP <- reactiveVal(NULL)
    result_G2_MF <- reactiveVal(NULL)
    result_G2_KEGG <- reactiveVal(NULL)

    observeEvent(geneListReactive(), {
      req(topGeneList(), downGeneList(), settingMAE())

      result_G1_CC(NULL)
      result_G1_BP(NULL)
      result_G1_MF(NULL)
      result_G2_CC(NULL)
      result_G2_BP(NULL)
      result_G2_MF(NULL)

      mae <- settingMAE()
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- c("G1", "G2")
      group1_fc_gene_profile <- topGeneList()
      group2_fc_gene_profile <- downGeneList()

      for (n in seq_along(groups_list)) {
        col <- groups_list[n]
        gene_list <- get(c("group1_fc_gene_profile", "group2_fc_gene_profile")[n])
        for (mode in c("CC", "BP", "MF")) {
          future_promise({
            start_time <- Sys.time()
            result <- go_enrich_dotplot(
              gene_list_ = unique(gene_list),
              save_path_ = NULL,
              save_filename_ = NULL,
              mode_ = mode,
              showCategory_ = 10
            )
            end_time <- Sys.time()
            list(
              r = result, c = col, m = mode,
              start_time = start_time,
              end_time = end_time,
              elapsed = as.numeric(difftime(end_time, start_time, units = "secs"))
            )
          }) %...>% {
            if (!is.null(.)) {
              var <- paste0(.$c, "_", .$m)
              cat(var, " started at:", as.character(.$start_time), "\n")
              cat(var, " ended at:", as.character(.$end_time), "\n")
              cat(var, " elapsed:", as.character(.$elapsed), " seconds\n")


              if (var == "G1_CC") result_G1_CC(.$r)
              if (var == "G1_BP") result_G1_BP(.$r)
              if (var == "G1_MF") result_G1_MF(.$r)
              if (var == "G2_CC") result_G2_CC(.$r)
              if (var == "G2_BP") result_G2_BP(.$r)
              if (var == "G2_MF") result_G2_MF(.$r)
            }
          }
        }
      }
    })


    observeEvent(geneListReactive(), {
      req(topGeneList(), downGeneList(), settingMAE())

      result_G1_KEGG(NULL)
      result_G2_KEGG(NULL)

      mae <- settingMAE()
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- c("G1", "G2")
      group1_fc_gene_profile <- topGeneList()
      group2_fc_gene_profile <- downGeneList()

      for (n in seq_along(groups_list)) {
        col <- groups_list[n]
        gene_list <- get(c("group1_fc_gene_profile", "group2_fc_gene_profile")[n])
        future_promise({
          start_time <- Sys.time()
          result <- kegg_enrich_dotplot(
            gene_list_ = unique(gene_list),
            save_path_ = NULL,
            save_filename_ = NULL,
            showCategory_ = 10
          )
          end_time <- Sys.time()
          list(
            r = result, c = col,
            start_time = start_time,
            end_time = end_time,
            elapsed = as.numeric(difftime(end_time, start_time, units = "secs"))
          )
        }) %...>% {
          if (!is.null(.)) {
            var <- paste0(.$c, "_", "KEGG")
            cat(var, " started at:", as.character(.$start_time), "\n")
            cat(var, " ended at:", as.character(.$end_time), "\n")
            cat(var, " elapsed:", as.character(.$elapsed), " seconds\n")

            if (var == "G1_KEGG") result_G1_KEGG(.$r)
            if (var == "G2_KEGG") result_G2_KEGG(.$r)
          }
        }
      }
    })


    output$G1_CC <- renderPlot({
      validate(need(result_G1_CC(), "Please update the DEG list, then wait while it loads..."))
      result_G1_CC()
    })

    output$G1_BP <- renderPlot({
      validate(need(result_G1_BP(), "Please update the DEG list, then wait while it loads..."))
      result_G1_BP()
    })
    output$G1_MF <- renderPlot({
      validate(need(result_G1_MF(), "Please update the DEG list, then wait while it loads..."))
      result_G1_MF()
    })
    output$G2_CC <- renderPlot({
      validate(need(result_G2_CC(), "Please update the DEG list, then wait while it loads..."))
      result_G2_CC()
    })
    output$G2_BP <- renderPlot({
      validate(need(result_G2_BP(), "Please update the DEG list, then wait while it loads..."))
      result_G2_BP()
    })
    output$G2_MF <- renderPlot({
      validate(need(result_G2_MF(), "Please update the DEG list, then wait while it loads..."))
      result_G2_MF()
    })
    output$G1_KEGG <- renderPlot({
      validate(need(result_G1_KEGG(), "Please update the DEG list, then wait while it loads..."))
      result_G1_KEGG()
    })
    output$G2_KEGG <- renderPlot({
      validate(need(result_G2_KEGG(), "Please update the DEG list, then wait while it loads..."))
      result_G2_KEGG()
    })

    outputOptions(output, "G1_CC", suspendWhenHidden = FALSE)
    outputOptions(output, "G1_BP", suspendWhenHidden = FALSE)
    outputOptions(output, "G1_MF", suspendWhenHidden = FALSE)
    outputOptions(output, "G2_CC", suspendWhenHidden = FALSE)
    outputOptions(output, "G2_BP", suspendWhenHidden = FALSE)
    outputOptions(output, "G2_MF", suspendWhenHidden = FALSE)
    outputOptions(output, "G1_KEGG", suspendWhenHidden = FALSE)
    outputOptions(output, "G2_KEGG", suspendWhenHidden = FALSE)

    observe({
      req(result_G1_CC(), result_G1_BP(), result_G1_MF(), result_G2_CC(), result_G2_BP(), result_G2_MF(), result_G1_KEGG(), result_G2_KEGG())
      shinyjs::enable("run_DEG")
      shinyjs::enable("dbBrowser1-selected_db")
    })

    observeEvent(geneListReactive(), {
      req(DEG_table())
      gseaFCModuleServer("gsea_up", DEG_table = DEG_table, direction = "up")
      gseaFCModuleServer("gsea_down", DEG_table = DEG_table, direction = "down")
    })
    observe({
      req(DEG_table(), wide_data(), maeColData())
      print(head(maeColData()))
      print(head(wide_data()))
      pcaModuleServer("pca1", wide_data(), maeColData())
    })
  }

  for_run <- shinyApp(ui = ui, server = server)
  runapp <- runApp(for_run)

  return(runapp)
}
